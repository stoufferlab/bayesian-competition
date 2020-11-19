require(tidyverse)
require(mvtnorm)
require(cubature)
require(Gmedian)


#This is saavedras function, it uses an interaction matrix to integrate using a multivariate normal distribution
#They assume growth rates have no constraints and an radius (R=1), but in their case this result is valid fro any value of R if you use a LV model
Omega_SA <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol =  1.88129e-25 )
  #calculate the prob density
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  return(d[1])
 
}


#This is our function. We do care about constraints on the growht rates and on the abundances.
#Given a growth rate of species "i" inside a Radius (R), and 
#constraints on the growthrates (rconstraints)
#and abundances (Nupper)
#It determines the value of rj
#It determines if the vector of ri,rj is feasible given the above. 
feasibility <- function(ri,alpha,R,branch=c("plus","minus"),rconstraints=NULL,Nupper=NULL){
  # growth rates are on a circle of radius R so given one we can determine the other
  rj <- switch(branch,
               plus  = sqrt(R^2 - ri^2),
               minus = -sqrt(R^2 - ri^2)
  )
  r <- c(ri,rj)
  
  # solve for the equilibrium given the interactions and the growth rate vector
  N <- solve(alpha) %*% r
  
  # check if N corresponds to feasibile equilibrium
  N_feasible <- (N > 0)
  
  # check whether N are within constraints
  if(!is.null(Nupper)){
    N_good <- (N <= Nupper)
  }else{
    N_good <- rep(TRUE,length(N))
  }
  
  # check whether growth rates are within constraints
  if(!is.null(rconstraints)){
    r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
  }else{
    r_good <- rep(TRUE,length(r))
  }
  
  return(prod(N_feasible, N_good, r_good))
}

#same function as above but vectorized. This doesnt really matter, I was just following the vignette and it said this was faster, but it gives the same result
vectorized_feasibility <- function(x, alpha, R, branch=c("plus","minus"),rconstraints=NULL,Nupper=NULL){
 
  matrix(sapply(
    x,
    function(x,  alpha, R, branch=c("plus","minus"),rconstraints=NULL,Nupper=NULL){
      rj <- switch(branch,
                   plus  = sqrt(R^2 - x^2),
                   minus = -sqrt(R^2 - x^2)
      )
      r <- c(x,rj)
      
      # solve for the equilibrium given the interactions and the growth rate vector
      N <- solve(alpha) %*% r
      
      # check if N corresponds to feasibile equilibrium
      N_feasible <- (N > 0)
      
      # check whether N are within constraints
      if(!is.null(Nupper)){
        N_good <- (N <= Nupper)
      }else{
        N_good <- rep(TRUE,length(N))
      }
      
      # check whether growth rates are within constraints
      if(!is.null(rconstraints)){
        r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
      }else{
        r_good <- rep(TRUE,length(r))
      }
      
      return(prod(N_feasible, N_good, r_good)) 
    }, alpha=alpha, R =R, branch=branch, rconstraints=rconstraints, Nupper=Nupper
  ) , ncol = ncol(x))
   
}


# A function to integrate over the whole sphere, using a determined R
omega <- function(R,alpha,rconstraints=NULL,Nupper=NULL){
  #nevermind how I call the omegas, it basically integrates over the 4 quadrants
  # calculate the "plus" branch of the integral for the alpha matrix
  omega_right <- adaptIntegrate(
    feasibility,
    lower=c(0),
    upper=c(R),
    alpha=alpha,
    R=R,
    branch="plus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  
  omega_left <- adaptIntegrate(
    feasibility,
    lower=c(-R),
    upper=c(0),
    alpha=alpha,
    R=R,
    branch="plus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  
  # calculate the "minus" branch of the integral for the alpha matrix
  omega_up <- adaptIntegrate(
    feasibility,
    lower=c(0),
    upper=c(R),
    alpha=alpha,
    R=R,
    branch="minus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  
  
  omega_down <- adaptIntegrate(
    feasibility,
    lower=c(-R),
    upper=c(0),
    alpha=alpha,
    R=R,
    branch="minus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  # the total integral is the sum of the four
  omega_alpha <- omega_right$integral + omega_left$integral + omega_down$integral + omega_up$integral
  

  # because of the limits of the integral on r, we can normalize directly
  omega <- omega_alpha / (4*R)
  
  return(omega)
}

#A function to integrate again, using the vectorized feasibility. It uses different integration methods but it gives the same result. again just for following the instructions on the vignette
vectorized_omega<-function(R,alpha,rconstraints=NULL,Nupper=NULL){
  
  omega_plus  <- cubintegrate(f=vectorized_feasibility,
                              lower = c(-R),
                              upper = c(R),
                              method = "cuhre",
                              alpha=alpha,
                              R=R,
                              branch="plus",
                              rconstraints=rconstraints,
                              Nupper=Nupper,
                              nVec = 1024)
  
  omega_minus  <- cubintegrate(f=vectorized_feasibility,
                               lower = c(-R),
                               upper = c(R),
                               method = "cuhre",
                               alpha=alpha,
                               R=R,
                               branch="minus",
                               rconstraints=rconstraints,
                               Nupper=Nupper,
                               nVec = 1024)
  # the total integral is the sum of the two
  omega_alpha <- omega_plus$integral + omega_minus$integral
  
  # because of the limits of the integral on r, we can normalize directly
  omega <- omega_alpha / (4*R)
  
  return(omega)
}

#Function to see how different is our integration from Saavedras.  
#We create an 2 by 2 alpha matrix with a determined values of intraspecific competition (diagonal)
#We see how the calculation of omegas changes as interspecific competition (off diagonal) increases and plot it
multiple_omega<-function(diagonal,off_diagonal, R, rconstraints=NULL, Nupper=NULL){
  
  alpha <- alpha<- diag(diagonal, ncol = 2, nrow = 2)
  inter_alphas <- seq(0, off_diagonal, 0.01)
  
  omegas <- sapply(inter_alphas, function(x, R, rconstraints, Nupper){
    
    alpha[1,2] <- x
    alpha[2,1] <- x
   
    
   vector_integration <- vectorized_omega(R=R,
                     alpha = alpha,
                     rconstraints = rconstraints, 
                     Nupper = Nupper )
   

   feasibility_integration <-  omega(R=R,
                                     alpha = alpha,
                                     rconstraints = rconstraints, 
                                     Nupper = Nupper )
  
   
   saavedra_integration <- Omega_SA(alpha)
   
   
    results <-  c( x,vector_integration, feasibility_integration, saavedra_integration) %>% as.data.frame()
   
    return( results)
    
  },R =R, rconstraints =rconstraints, Nupper=Nupper   )
  
   omegas_together <- do.call(rbind, omegas) %>% as.data.frame()
   
   colnames(omegas_together)<- c("inter_alpha", "vec_omega", "omega", "saavedra")
    
   
   plot(omegas_together$inter_alpha, 
        omegas_together$vec_omega,
        type = "l",
        col = "dodgerblue",
        lwd = 2,
        xlab = "Intraspecific alpha",
        ylab = "Feasibility Domain")
   
   lines(omegas_together$inter_alpha,
         omegas_together$omega,
         col = "seagreen",
         lwd = 2)
   
   
   lines(omegas_together$inter_alpha,
         omegas_together$saavedra,
         col = "firebrick",
         lwd = 2)
   return(omegas_together)
  
}



################################################################################# 
# If there is no constraints on the growth rates or abundances, our results SHOULD be the same as Saavedras. But they are not! check this out

rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf)
)
Nupper <- c(
  i = Inf,
  j = Inf
)
#It doesnt really matter if we change R because we have no constraints but just for the sake of it
R<-1

#If there is no intraspecific competition (or very low, to avoid having a singlular matrix) the feasibility domain does not change as intraspecific competition grows. Our calculations are equal to Saavedras
multiple_omega( diagonal = 0.0001,
  off_diagonal = .9,
               R = 1, 
               rconstraints = rconstraints, 
               Nupper =  Nupper)

#If intraspecific competition is = 1 ...
#Red lines are Saavedras
#Green and blue lines are our integration
multiple_omega( diagonal = 1,
                off_diagonal = .9,
                R = 1, 
                rconstraints = rconstraints, 
                Nupper =  Nupper)


multiple_omega( diagonal = .513,
                off_diagonal = .9,
                R = 1, 
                rconstraints = rconstraints, 
                Nupper =  Nupper)


library(cubature)
library(Gmedian)

#Checks if a particular combination of growth rates, inside  a Radius, is feasible with the interaction matrix and constraints
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

#perform integration using the function above to determine the size of the feasibility domain
#For a particular R, it performs the integration of feasibility over its volume, returns normalized omega
omega <- function(R,alpha,rconstraints=NULL,Nupper=NULL){
  # calculate the "plus" branch of the integral for the alpha matrix
  omega_plus <- adaptIntegrate(
    feasibility,
    lower=c(-R),
    upper=c(R),
    alpha=alpha,
    R=R,
    branch="plus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  # calculate the "minus" branch of the integral for the alpha matrix
  omega_minus <- adaptIntegrate(
    feasibility,
    lower=c(-R),
    upper=c(R),
    alpha=alpha,
    R=R,
    branch="minus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  # the total integral is the sum of the two
  omega_alpha <- omega_plus$integral + omega_minus$integral
  
  # because of the limits of the integral on r, we can normalize directly
  omega <- omega_alpha / (4*R)
  
  return(omega)
}

#Calculate Omega for different values of R
omega_radius<-function(R_vals, alpha,rconstraints=NULL, Nupper=NULL){

  
  # calculate the proportional feasibility domain across the different lengths of R
  omega_vals <- sapply(
    R_vals,
    omega,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  
  plot(R_vals, omega_vals, type='l', xlab=expression(italic(R)), ylab=expression(Omega))
  return(data.frame("R"=R_vals, "Omegas"= omega_vals))
}

# sample values of R that are feasible
r_feasible<-function(R_vals,alpha, rconstraints=NULL,Nupper=NULL){

  r_sample <- t(sapply(
    seq_len(1000),
    function(x,R_vals,alpha,rconstraints,Nupper){
      while(TRUE){
        #you sample a value of R
        R <- sample(R_vals,1)

        theta <- runif(1,0,2*pi)

        ri <- R*cos(theta)
        branch <- ifelse(sin(theta)>0,"plus","minus")
        feas <- feasibility(ri,alpha,R,branch,rconstraints,Nupper)
        if(feas){
          return(c(R=R,theta=theta,ri=ri,rj=R*sin(theta)))
        }
      }
    },
    R_vals=R_vals,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  ))
  plot(0,0,
       xlim=c(-range(R_vals)[2],range(R_vals)[2]),
       ylim=c(-range(R_vals)[2],range(R_vals)[2]),
       type='n',
       xlab=expression(italic(r[i])),
       ylab=expression(italic(r[j]))
  )
  abline(h=0,lty='dashed',lwd=1.5)
  abline(v=0,lty='dashed',lwd=1.5)
  apply(r_sample, MARGIN=1,
        function(x) {
          segments(0,0,x["ri"],x["rj"],col=grey(0.75))
        }
  )
  

  # compute and plot the median vector (which appears to end up a bit short for some reason)
  medianR <- Gmedian(r_sample[,c("ri","rj")])
  segments(0,0,medianR[1],medianR[2],lwd=2,col="aquamarine4")
  
  return(r_sample)
}

feasibility_wrapper<-function(max_R, alpha, rconstraints=NULL, Nupper=NULL){
  R_vals <- seq(0.0001, max_R, length.out = 100)
  #we calculate omega values for different values of R
  #R_vals <- c(0.0001,1:100)
  omega_values <- omega_radius(R_vals = R_vals,
                               alpha = alpha,
                               rconstraints = rconstraints,
                               Nupper = Nupper)
  #We sample the values of R that are feasible
  growth_values<- r_feasible(R_vals = R_vals,
                             alpha = alpha,
                             rconstraints = rconstraints,
                             Nupper = Nupper  )
  #We average over omega_values to get the mean feasibility domain
  mean_omega<- apply(omega_values,2,mean)[2]
  #We get the median of the growth vectors to calculate the centroid
  
  r_centroid<- Gmedian(growth_values[,c("ri","rj")])
  #return(growth_values)
  return(list("Omega"=mean_omega, "r_centroid"=r_centroid))
  
}


theta <- function(r_c,r){
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}


alpha <- diag(2)
alpha[1,2] <- 0.80
alpha[2,1] <- 0.223
# are there values of growth rate that we do not wish to consider?
rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c( Inf, Inf)
)

# are there upper bounds on N that we do not consider biologically realistic?
Nupper <- c(
  i = Inf,
  j = Inf
)

feasibility_wrapper(max_R = 100,
                    alpha = alpha,
                    rconstraints = rconstraints,
                    Nupper = Nupper)

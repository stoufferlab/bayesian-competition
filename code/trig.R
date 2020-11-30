require(cubature)
require(tidyverse)



#previous functions used by saaveda et al. 
Omega_SA <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha )
  #calculate the prob density
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  return(d[1])
  
}
omega_analytical <- function(alpha){
  circle <- 2/pi 
  num <- (alpha[1,1] * alpha[2,2]) - (alpha[1,2] * alpha[2,1])
  sp1 <- sqrt(alpha[1,1]^2 + alpha[2,1]^2)
  sp2 <- sqrt(alpha[1,2]^2 + alpha[2,2]^2)
  denom <- sp1 * sp2
  domain <- asin((num/denom))
  feas <- circle*domain
  #return(feas)
  return(feas/4)
  
}
calculate_area<-function(R=1, alpha){
  
  #test<- r_feasible(R=R, alpha = alpha, make_plot =TRUE)
  
  
  #the slopes
  m1 <- alpha[2,2]/ alpha[1,2]
  m2 <- alpha[2,1]/ alpha[1,1]
  
  #The inclination angle
  theta1 <- atan(m1)
  theta2 <- atan(m2)
  
  #their difference
  final_theta <- abs(theta1 - theta2)
   
  #The final angle in radians
  angle <- (2 *pi) - final_theta
  
  #we convert it to degrees
  angle_degrees <- (angle *180) / pi
  
  
  #DrawArc(x=0,y=0, rx=R, theta.1 = theta1, theta.2=theta2)
  
  #The area of the full circle 
  area_circle<- pi * (R^2)
  
  #The area of the sector outside the feasibility domain
  area_sector <- (area_circle * angle_degrees ) / 360
  
  #The area of the feasibility domain
  feasibility_area <-area_circle- area_sector
 
 
 ## but it is not the final area we want, we want the proportion of feasible space
 # the maximum space is the area of the circle
 
 proportion_area <- feasibility_area / area_circle
 
 
 return(proportion_area)
  
}


#given a value of theta (in radians) it calculates the coresponding ri and rj, and checks if they are feasible
feasibility_theta <- function(theta,alpha, R,rconstraints=NULL,Nupper=NULL){
  
  
  #The magnitude of the vector is its R
  # Theta can go from 0 to 2pi
  
  ri <- R * cos(theta)
  rj <- R * sin(theta)
  
  r <- c(ri,rj)
  # solve for the equilibrium given the interactions and the growth rate vector
  N <- solve(alpha) %*% r
  
  # check if N corresponds to feasibile equilibrium
  N_feasible <- (N > 0)
  
  if(!is.null(Nupper)){
    N_good <- (N <= Nupper)
  }else{
    N_good <- rep(TRUE,length(N))
  }

  # check whether growth rates are within constraints
  if(!is.null(rconstraints)){
    r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
    #print(r_good)
  }else{
    r_good <- rep(TRUE,length(r))
  }

  feasible <- prod(N_feasible* N_good* r_good)
  
   # if(feasible){
   #   segments(0,0, ri, rj, col="dodgerblue")
   # }

  return(feasible)
 
}

#a numerical integration for all the values of theta in a circle of radius R
numerical_theta<- function(alpha, R, rconstraints=NULL, Nupper=NULL){
  thetas <- seq(0 , 2*pi, .01)
  
  numerical<- sapply(thetas, function(x, alpha, R, rconstraints, Nupper){
    pi <-feasibility_theta(theta = x,
                           alpha = alpha,
                           R=R,
                           rconstraints = rconstraints,
                           Nupper= Nupper)
  }, alpha = alpha, R=R, rconstraints=rconstraints, Nupper=Nupper)
  
  
  integration <-  sum(numerical)/length(thetas)
  return(integration)
  
}

# We integrate from 0 to the maximum value of R 
integrate_radii <- function( alpha, R_max, rconstraints=NULL, Nupper=NULL){
 
   all_radii <- cubintegrate(f=numerical_theta, 
                      lower = 0, 
                      upper = R_max,
                      method = "hcubature",
                      alpha = alpha,
                      rconstraints=rconstraints,
                      Nupper=Nupper)
  
  return(all_radii$integral)
} 

#check out how our integration changes with competition coeff
multiple_omega<-function(diagonal,off_diagonal, R_max, rconstraints=NULL, Nupper=NULL){
  
  alpha <- alpha<- diag(diagonal, ncol = 2, nrow = 2)
  inter_alphas <- seq(0, off_diagonal, 0.01)
  
  omegas <- sapply(inter_alphas, function(x, R_max, rconstraints, Nupper){
    
    alpha[1, 2] <- x
    alpha[2, 1] <- x
    
    saavedra_integration <- Omega_SA(alpha)
    analytical_intecration <- omega_analytical(alpha)
    numerical_integration <- integrate_radii(
      alpha = alpha,
      R_max = R_max,
      rconstraints = rconstraints,
      Nupper = Nupper )
    
    
    results <-  c( x,saavedra_integration, analytical_intecration, numerical_integration) %>% as.data.frame()
    
    return( results)
    
  },R_max=R_max, rconstraints =rconstraints, Nupper=Nupper   )
  
  omegas_together <- do.call(rbind, omegas) %>% as.data.frame()
  
  colnames(omegas_together)<- c("inter_alpha" ,"saavedra", "analytical","numerical")
  
  
  plot(omegas_together$inter_alpha,
       omegas_together$saavedra,
       col = "firebrick",
       pch=18)
  
  points(omegas_together$inter_alpha,
         omegas_together$numerical,
         col= "dodgerblue",
         pch=16)
  
  lines(omegas_together$inter_alpha,
        omegas_together$analytical,
        lwd=1.5,
        col="darkgoldenrod")
  
  return(omegas_together)
  
}



rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf)
)
Nupper <- c(
  i = Inf,
  j = Inf
)


test<-multiple_omega( diagonal = 1,
                    off_diagonal = .9,
                    R = 10, 
                    rconstraints = rconstraints, 
                    Nupper =  Nupper)



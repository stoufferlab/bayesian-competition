
require(tidyverse)
require(mvtnorm)
require(cubature)


#previous functions used by saaveda et al. These  functions give the exact same result for a given alpha and R.
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
# And the characteristics of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}
theta <- function(r_c,r){
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}
test_feasibility_saavedra <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}

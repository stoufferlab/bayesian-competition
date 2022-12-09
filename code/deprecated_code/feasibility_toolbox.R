require(tidyverse)
require(mvtnorm)

 Omega <- function(alpha){
   n <- nrow(alpha)
   Sigma <-solve(t(alpha) %*% alpha, tol =  1.88129e-25 )
   #calculate the prob density
   d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
   #out <- log10(d[1]) + n * log10(2)
  # return( d[1]^(1 / n))
  return(d[1])
 }


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


test_feasibility <- function(alpha,r){
  out <- all(solve(alpha,r)>0)
  return(out)
}
Omega_with_inequality_constraints <- function(alpha, inequality){
  beta <- alpha
  for(j in 1:ncol(alpha)){
    beta[,j] <- -alpha[,j] / (-alpha[,j] %*% inequality)
  }
  abs(det(beta))
}



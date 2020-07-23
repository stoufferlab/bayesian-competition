
require(mvtnorm)

# Omega <- function(alpha){
#   n <- nrow(alpha)
#   Sigma <-solve(t(alpha) %*% alpha)
#   d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
#   out <- log10(d[1]) + n * log10(2)
#   return(d[1])
# }




r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}

theta <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}


# Omega<-function(alpha){
#   num<- (alpha[1,1]*alpha[2,2]) - (alpha[1,2]*alpha[2,1] )
#   denom<- sqrt( alpha[1,1]^2 + alpha[2,1]^2) * sqrt(alpha[1,2]^2 + alpha[2,2]^2 )
#   a<-asin(num/denom)
#   omega<-(2/pi)*a
#   return(omega)
# }



Omega <- function(alpha) {
  S <- nrow(alpha)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
  #   if (length(which(diag(alpha) == 0)) == 0) {
  #     Sigma <- chol2inv(alpha, size = NCOL(alpha), LINPACK = FALSE)
  #     return(omega(S, Sigma))
  #   }
  #   else {
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (f(alpha) == FALSE) {
    return(0)
  }
  else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
  #   }
  # }
}


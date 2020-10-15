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




calculate_constrained_domain <- function(alpha, constraints){
  Nsample <- 1e5
  no_constraints<- list(c(-1e5,1e5),c(-1e5,1e5))
  
  # We sample N number of growth rates when there is no constraints
  r <- no_constraints %>% 
    map(~runif(Nsample, min = .[1], max = .[2])) %>% 
    unlist() %>% 
    matrix(ncol = Nsample, byrow = TRUE)
  
  #We test and save which ones are feasible
  feasibility <- c()
  r1 <- c()
  r2 <- c()
  for(i in 1:Nsample){
    
    feasibility[i] <-test_feasibility(alpha,r[,i])
    r1[i]<-r[,i][1]
    r2[i]<-r[,i][2]
    
  }
  
  full_domain<-cbind(feasibility,r1,r2) %>% as.data.frame()
  feasibility_domain <- filter(full_domain, feasibility ==1)
  
  #Now we sample the feasibility domain for our constraints
  feasibility_domain_constrained<- filter(feasibility_domain, constraints[1] & constraints[2] )
  
  omega_prime<-sum(feasibility_domain_constrained$feasibility) / Nsample
  r_centroid_prime<-apply(feasibility_domain_constrained,2,mean)[2:3]
  return(c(omega_prime, r_centroid_prime))
  
}

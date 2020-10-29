require(tidyverse)
#Consider that we estimate the parameters of both species with a B-H model
#we have the following interaction matrix
#Sp.1 facilitates species two but limits itself, sp.2 limits itself and species 1.
alpha<-matrix(c(1,0, 0, 1), nrow = 2, ncol = 2)
#We can calculate the feasibility domain by sampling growth rates, asuming no "constraints".
#Lets also assume that Infinity is a big number like 1e6

no_constraints<- list(c(-1e6,1e6),c(-1e6,1e6))
Nsample <- 1e5
# We sample N number of growth rates
r <- no_constraints %>% 
  map(~runif(Nsample, min = .[1], max = .[2])) %>% 
  unlist() %>% 
  matrix(ncol = Nsample, byrow = TRUE)


#And test and save which ones are feasible given our alpha
test_feasibility <- function(alpha,r){
  out <- all(solve(alpha,r)>0)
  return(out)
}

feasibility <- c()
r1 <- c()
r2 <- c()
for(i in 1:Nsample){

  feasibility[i] <-test_feasibility(alpha,r[,i])
  r1[i]<-r[,i][1]
  r2[i]<-r[,i][2]
  
}

#Originally you would calculate Omega by:
omega<- mean(feasibility)

#But we know that is an over estimation of the feasibility domain, because there are some constraints on the possible values of the growth rates.
#So we can sample inside the feasibility domain to see which proportion of it is compatible with our constraints. First we define the feasibility domaain

full_domain<-cbind(feasibility,r1,r2) %>% as.data.frame()
feasibility_domain <- filter(full_domain, feasibility ==1)

#Sp.1 and 2 can only have growth rates that go from -1, to Inf. So lets take out any combination that is outside those limits
feasibility_domain_constrained<- filter(feasibility_domain, r1 >= -1 & r2 >= -1 )

#And we can estimate now Omega prime as the new proportion of those values that fall inside the feasibility domain AND are within our constraints
omega_prime<-sum(feasibility_domain_constrained$feasibility) / Nsample


#This way we guarantee that Omega prime is the same size or smaller than Omega.














Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol =  1.88129e-25 )
  #calculate the prob density
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  #out <- log10(d[1]) + n * log10(2)
  # return( d[1]^(1 / n))
  return(d[1])
}






Omega_limits <- function(alpha,lower_sp_1, upper_sp_1, lower_sp_2, upper_sp_2 ){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol =  1.88129e-25 )
  #calculate the prob density
  r1_low<- lower_sp_1 / sum(alpha[,1]) 
  r1_up <- upper_sp_1 / sum(alpha[,1]) 
  
  r2_low<- lower_sp_2 / sum(alpha[,2]) 
  r2_up <- upper_sp_2 / sum(alpha[,2]) 
  
 
   d <- pmvnorm(lower = c(r1_low, r2_low), upper = c(r1_up,r2_up), mean = rep(0,n), sigma = Sigma)
  #out <- log10(d[1]) + n * log10(2)
  # return( d[1]^(1 / n))
  return(d)
}




Omega_2sp<-function(alpha){
  p <- 2/pi
  a  <- (alpha[1,1]*alpha[2,2]) - (alpha[1,2]*alpha[2,1])
  b  <- sqrt((alpha[1,1]^2) + (alpha[2,1]^2))
  c  <- sqrt((alpha[1,2]^2) + (alpha[2,2]^2))
  
  d <- b*c
  e <-a/d
  angle<- p*asin(e)
  return(angle)
}

r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}
projection_2sp<- function(alpha,x){
  
 
  
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[, 1]
  v2 <- alpha_n[, 2]
  vc <- (v1 + v2)
  vc <- vc / sqrt(sum(vc ^ 2))
  
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  #v3 <- v3/sum(v3)
  vc <- vc / sum(vc)
  
  
  
  
  plot(
    0,
    axes = F,
    col = 'grey50',
    type = 'l',
    lwd = 2,
    xlim = c(-x,x),
    ylim = c(-x,x),
    xlab = "",
    ylab = "",
    cex.lab = 1.5,
    cex.main = 2,
    adj  = 0
  )
  lines(c(0, 0), c(0, x), col = 'grey50', lwd = 2)
  lines(c(0, 0), c(0, -x), col = 'grey50', lwd = 2)
  lines(c(-x,0), c(0, 0), col = 'grey50', lwd = 2)
  lines(c(0,x), c(0, 0), col = 'grey50', lwd = 2)
  #####
  
  
  
  vcP <- v1 / sqrt(sum(v1 ^ 2)) + v2 / sqrt(sum(v2 ^ 2))
  
  vcP <- vcP / sum(vcP)
  
  
  points( v1[1], v1[2], col='mediumseagreen',pch=16,cex=1)
  points( v2[1], v2[2], col='mediumseagreen',pch=16,cex=1)
  points( vcP[1], vcP[2], col="orange1", pch=16,cex=1)
  
  
  
  # 
  # upper <- alpha[2, 2] / alpha[1, 2]
  # lower <- alpha[2, 1] / alpha[1, 1]
  # 
  # xx <- seq(0, 2, .01)
  # yy <- seq(0, 2, .01)
  # lines(xx,
  #       xx * lower,
  #       col = 'mediumseagreen',
  #       lty = 1,
  #       lwd = 2)
  # lines(yy,
  #       yy * upper,
  #       col = 'mediumseagreen',
  #       lty = 1,
  #       lwd = 2)
  # 
  
  
  #growth <- r2 / r1
  #zz <- seq(0, vcP[1]-.01, 0.01)
  # lines(zz,
  #       zz * growth,
  #       col = "tomato3",
  #       lwd = 1.5,
  #       lty = 2)
  #arrows(0,0, max(zz), max(zz*growth), col="tomato3", lwd=2)
  
  #print(max(zz*growth))
 # centroid <- r_centroid(alpha)
  #points(centroid[1, ], centroid[2, ], col = "orange1", pch = 18)
  
  #slope_centroid <- (centroid[2, ] - vcP[2]) / (centroid[1, ] - vcP[1])
  
  #lines(zz, slope_centroid * zz, lwd = 2, col = "orange1")
  
  
  
}

alpha<-matrix(c(1,-0.7, -0.7, 1), nrow = 2, ncol = 2)


Omega(alpha_m)
Omega_2sp(alpha_m)
projection_2sp(alpha_m)


Omega_limits(alpha = alpha_m,lower_sp_1 = 0,upper_sp_1 = Inf,lower_sp_2 = 0,upper_sp_2 = Inf)

lines(c(-1,-1),c(-1,4), lty="dashed")
lines(c(-1,4),c(-1,-1), lty="dashed")

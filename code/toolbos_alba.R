
require(mvtnorm)

Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(exp(out)) 
  
}
Omega_song <- function(alpha) {
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

projection_2sp<- function(alpha,r1,r2,titulo){
  
  #graphics.off()
  #pdf(file=file_name, width = 8, height = 8, useDingbat=FALSE, family = "ArialMT")
 # par(mar = c(3,3,3,3))
  
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
 # v3 <- alpha_n[,3]
 # vc <- (v1 + v2 + v3)
  vc <- (v1 + v2 )
  vc <- vc / sqrt(sum(vc^2))
  
  v1 <- v1/sum(v1)
  v2 <- v2/sum(v2)
  #v3 <- v3/sum(v3)
  vc <- vc/sum(vc)
  

  

  plot(c(0,2),c(0,0),axes=F,col='grey50',type='l',lwd=2,xlim=c(0,1),ylim=c(0,1), main=titulo, xlab="Intrinsic growth rate vero (r1)",ylab="Intrinisc growth rate trcy (r2)")
  lines(c(0,0), c(0,2), col='grey50', lwd=2)
  #####
  
  
  
  vcP <- v1/sqrt(sum(v1^2)) + v2/sqrt(sum(v2^2))

  vcP <- vcP / sum(vcP)
  
  
  points( v1[1], v1[2], col='dodgerblue',pch=16,cex=1.5)
  points( v2[1], v2[2], col='dodgerblue',pch=16,cex=1.5)
  points( vcP[1], vcP[2], col="orange3", pch=16,cex=1.5 )
  
   upper<-alpha[2,2]/alpha[1,2]
   lower<-alpha[2,1]/alpha[1,1]
  
   xx<-seq(0,v1[1],.01)
   yy<-seq(0,v2[1],.01)
   lines(xx,xx*lower,col='mediumseagreen',lty=2,lwd=2)
   lines(yy,yy*upper,col='mediumseagreen',lty=2,lwd=2)
   
   growth<-r2/r1
   zz<-seq(0,vcP[1],0.01)
   lines(zz,zz*growth, col="brown",lwd=2)
   
   centroid<-r_centroid(alpha)
   points(centroid[1,],centroid[2,],col="orange1",pch=18)
  
}

projection_posterior<-function(alpha,r1,r2,col_alpha){
  omega<-Omega(alpha)
  feasibility<-test_feasibility(alpha,c(r1,r2))
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  vc <- (v1 + v2 )
  vc <- vc / sqrt(sum(vc^2))
  v1 <- v1/sum(v1)
  v2 <- v2/sum(v2)
  vc <- vc/sum(vc)
  vcP <- v1/sqrt(sum(v1^2)) + v2/sqrt(sum(v2^2))
  vcP <- vcP / sum(vcP)
  
  cblue<-rethinking::col.alpha("dodgerblue",alpha = col_alpha)
  corange<-rethinking::col.alpha("orange3",alpha = col_alpha)
  cgreen<-rethinking::col.alpha("mediumseagreen",alpha = col_alpha)
  cbrown<-rethinking::col.alpha("brown",alpha = col_alpha)
   
   points( v1[1], v1[2], col=cblue,pch=16,cex=1.5)
   points( v2[1], v2[2], col=cblue,pch=16,cex=1.5)
   points( vcP[1], vcP[2], col=corange, pch=16,cex=1.5 )
   
   upper<-alpha[2,2]/alpha[1,2]
   lower<-alpha[2,1]/alpha[1,1]
   
   xx<-seq(0,v1[1],.01)
   yy<-seq(0,v2[1],.01)
   lines(xx,xx*lower,col=cgreen,lty=2,lwd=1)
   lines(yy,yy*upper,col=cgreen,lty=2,lwd=1)
  
  growth<-r2/r1
   zz<-seq(0,vcP[1],0.01)
   lines(zz,zz*growth, col=cbrown,lwd=2)
   
    return(c(omega,feasibility))
}





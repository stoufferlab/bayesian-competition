source('code/toolbox_coexistence.R')
library(plotrix)
#this works for two species only

 a <- matrix(c(10,5,5,-1),2,2)
 r <- c(100,10)
 (Nstar <- solve(a) %*% r)

BV_growth<-function(g,s,lambda){
  num<- g*lambda
  den<- 1- ((1-g)*s)
  r<-   (num/den) - 1
}
#a11, a21, a12, a22
#-, + , +, +
#+, -, + , +
#+, +, - , +
#+, +, +, -
#-, - , +, +
#-, + , -, +

#-, + , +, -
#

#+,-,  - , +
#+,- , + , -



alpha <- matrix(c(.9,.4,.5,1),2,2)
domain<-function(alpha){
  upper<-alpha[2,2] / alpha[1,2]
  lower<- alpha[2,1] / alpha[1,1]
  centroid<-r_centroid(alpha)
  niche<- Omega(alpha)
  n<-100
  xx<-c(0,0)
  yy<-c(-1,1)
  par(mar = c(0,0,0,0))
  plot(xx,yy,xlim=c(-1,1),ylim=c(-1,1),type="l",xlab="Intrinsic growth rate sp 1", ylab="intrinsic growth rate sp 2")
  lines(yy,xx)
  
  arrows(0,0,2,upper*2,lwd=2, col='mediumseagreen')
  arrows(0,0,2,lower*2, lwd=2, col='mediumseagreen')
  
  points(centroid[1], centroid[2],pch=20, col="darkorange")
  
}
domain(alpha)


alpha <- matrix(c(1,0.0001,0.0001,1),2,2)

two_sp_omega<-function(alpha){
  num<- (alpha[1,1]*alpha[2,2]) - (alpha[2,1]-alpha[1,2] )
  denom<- sqrt( alpha[1,1]^2 + alpha[2,1]^2) * sqrt(alpha[1,2]^2 + alpha[2,2]^2 )
  a<-asin(num/denom)
  omega<-(2/pi)*a
  omega
}

two_sp_omega(alpha)
Omega(alpha)


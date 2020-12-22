
require(tidyverse)
require(mvtnorm)
require(cubature)
require(Gmedian)

#previous functions used by saaveda et al. These 3 functions give the exact same result for a given alpha and R.
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
#But what happens when we have constraints?


inverse_matrix<-function(alpha){
   
  deteminant_alpha <- (alpha[1,1] * alpha[2,2]) - (alpha[2,1] * alpha[1,2])
   
 inverse_det  <- 1 / deteminant_alpha
  
  
  adjugate <- matrix( NA,nrow=2,ncol=2)                  
  adjugate[1,1] <- alpha[2,2]
  adjugate[2,2] <- alpha[1,1]
  
  adjugate[1,2] <- -alpha[1,2]
  adjugate[2,1] <- -alpha[2,1]
  
  
  ii <- inverse_det * adjugate
  
  return(ii)
  
  
}
#given a value (or multiple values) of theta (in radians) it calculates the coresponding ri and rj, and checks if they are feasible

feasibility_theta <- function(theta_seq,alpha, R,rconstraints=NULL,Nupper=NULL){
  #we do not use solve to save time
  
  #we integrate over values of theta
  results <- sapply(theta_seq, function(theta, R, rconstraints, Nupper){
    #print(theta)
    ri <- R * cos(theta)
    rj <- R * sin(theta)
    r <- c(ri, rj)
    
    #we check first if the growth rates are within our constraints
    if (!is.null(rconstraints)) {
      r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
      r_feasible <- all(r_good)
      
    }else{
      r_feasible <- c(TRUE)
    }
    
    #if they are not, we do not bother to solve for abundances
    if (!r_feasible) {
      return(0)
    } else{
      #solve fo abundances
      inverse_alpha <- inverse_matrix(alpha)
      N1 <- (inverse_alpha[1, 1] * ri) + (inverse_alpha[1, 2] * rj)
      N2 <- (inverse_alpha[2, 1] * ri) + (inverse_alpha[2, 2] * rj)
      #check their feasibility
      N <- c(N1, N2)
      N_feasible <- (N > 0)
      N_feasible <- all(N_feasible)
      
      if(!N_feasible){
        return(0)
      }else{
         if (!is.null(Nupper)) {
           #we check that abundances are below the max abundance
           N_good <- (N <= Nupper) 
           N_good <- all(N_good)
           
         } else{
           N_good <- TRUE
       }
         return(N_good)
       
      }
      
    }
  }, R=R, rconstraints = rconstraints, Nupper=Nupper)
  
  
  return(results)
  
  
}

#vectorized integration of theta . We do a manual integration because sometimes the area is too small and adaptive integration misses it!
integrate_theta <-function( R_seq,alpha,rconstraints=NULL,Nupper=NULL ) {
  #we break down a circle into very small parts
  thetas <- seq(0 , 2*pi, length.out = 1000)
  
  results <- matrix(sapply(R_seq, function(R, alpha, rconstraints, Nupper){
     
     area<- feasibility_theta(theta_seq = thetas,
                              alpha = alpha,
                              R = R,
                              rconstraints = rconstraints,
                              Nupper= Nupper)
 
    #print(sum(area)/length(thetas))
    return(sum(area)/length(thetas))
  }, alpha=alpha,
  rconstraints = rconstraints,
  Nupper = Nupper))
  
 # results<- as.matrix(results)
  
  return(results)
  
}

#We integrate over different values of R to check how constraints affect the feasibility domain as R increases
integrate_radii <- function(alpha, R ,rconstraints=NULL,Nupper=NULL){

   multiple_R <- hcubature(f=integrate_theta,
                      lowerLimit  = 0, 
                      upperLimit  = R,
                      vectorInterface = TRUE,
                      tol = 1e-3,
                      alpha = alpha,
                      rconstraints= rconstraints,
                      Nupper= Nupper)
   
   return(multiple_R$integral/R)
   
   


}




#Sample growth rates inside the feasibility domain to determine its median or the center of the feasibility domain
r_feasible<-function(alpha, rconstraints=NULL, Nupper=NULL,R_max ,make_plot=FALSE){
  R_vals <- seq(0, R_max, length.out = 2000)
  
  #We sample values of R that are feasible to calculate their median, or the area in the center
  r_sample <- t(sapply(
    seq_len(1000),
    function(x,R_vals,alpha,rconstraints,Nupper){
      while(TRUE){
        R <- sample(R_vals,1)
        theta <- runif(1,0,2*pi)
        feas <- feasibility_theta(theta = theta,
                                  R= R,
                                  alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper)
        if(feas){
          return(c(R=R,theta=theta,ri=R*cos(theta),rj=R*sin(theta)))
        }
      }
    },
    R_vals=R_vals,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  )) %>% as.data.frame()
  
  
  medianR <- Gmedian(r_sample[,c("ri","rj")])

  if(make_plot){
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
  }
  return(medianR)
 # return(r_sample)
}

#Check if ourparticular combination of growth rates is feasible
check_feasibility <- function(r,alpha,rconstraints=NULL,Nupper=NULL){
  inverse_alpha <- inverse_matrix(alpha)
 
  #we check first if the growth rates are within our constraints
    if (!is.null(rconstraints)) {
      r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
      r_feasible <- all(r_good)
       }else{
      r_feasible <- c(TRUE,TRUE)
    }
    
    #if they are not, we do not bother to solve for abundances
    if (!r_feasible) {
      return(0)
    }else{
      #solve fo abundances
      N1 <- (inverse_alpha[1, 1] * r[1]) + (inverse_alpha[1, 2] * r[2])
      N2 <- (inverse_alpha[2, 1] * r[1]) + (inverse_alpha[2, 2] * r[2])
      #check their feasibility
      N <- c(N1, N2)
      N_feasible <- (N > 0)
      N_feasible <- all(N_feasible)
      
      if(!N_feasible){
        return(0)
      }else{
        if (!is.null(Nupper)) {
          N_good <- (N <= Nupper) 
          N_good <- all(N_good)
        } else{
          N_good <- TRUE
        }
        return(N_good)
      }
      
    }
  
}


#How far away are the observed growth rates from the center of the feasiblity domain?
calculate_distance<-function(center, r){
  dist<- sqrt((r[1]- center[1])^2 + (r[2]- center[2])^2)
  return(dist)
}




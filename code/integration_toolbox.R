
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
feasibility_theta <- function(theta_seq,alpha, R,rconstraints=NULL,Nupper=NULL){
  
  
  results <- sapply(theta_seq, function(theta, alpha, R, rconstraints, Nupper){
    #print(theta)
    ri <- R * cos(theta)
    rj <- R * sin(theta)
    r <- c(ri, rj)
    
    #we check first if the growth rates are within our constraints
    if (!is.null(rconstraints)) {
      r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
      r_good <- prod(r_good)
      
    } else{
      r_good <- TRUE
    }
    
    #if they are not, we do not bother to solve for abundances
    if (!r_good) {
      #  segments(0,0,ri,rj, col="firebrick")
      return(0)
    } else{
      # solve for the equilibrium given the interactions and the growth rate vector
      N <- solve(alpha) %*% r
      # check if N corresponds to feasibile equilibrium
      N_feasible <- (N > 0) %>% prod()
      
      #check if they are within our bounds of abundances
      if (!is.null(Nupper)) {
        N_good <- (N <= Nupper) %>% prod()
      } else{
        N_good <- TRUE
      }
      
      feasible <- prod(N_feasible * N_good)
      
      # if(feasible){
      #   segments(0,0,ri,rj, col="dodgerblue")
      # }
      
      return(feasible)
      
    }
  },alpha=alpha, R=R, rconstraints = rconstraints, Nupper=Nupper)
  
  
  return(results)
  
  
}


#vectorized integration of theta 
integrate_theta <-function( R_seq,alpha,rconstraints=NULL,Nupper=NULL ) {
  thetas <- seq(0 , 2*pi, length.out = 1000)
  
  results <- sapply(R_seq, function(R, alpha, rconstraints, Nupper){
     
     area<- feasibility_theta(theta_seq = thetas,
                              alpha = alpha,
                              R = R,
                              rconstraints = rconstraints,
                              Nupper= Nupper)
 
     # area<-integrate(f=feasibility_theta,
     #                lower = 0, 
     #                upper = 2*pi,
     #                alpha = alpha,
     #                R = R,
     #                rconstraints = rconstraints,
     #                Nupper= Nupper)
     # 
   
    return(sum(area)/length(thetas))
  }, alpha=alpha,
  rconstraints = rconstraints,
  Nupper = Nupper)
  return(results)
  
}


integrate_radii <- function(alpha, R ,rconstraints=NULL,Nupper=NULL){
  
   multiple_R <- integrate( f= integrate_theta,
                            lower = 0,
                            upper = R,
                            alpha = alpha,
                            rconstraints= rconstraints,
                            Nupper= Nupper)
   return( multiple_R[1]$value/R)

  # 
  # R_seq <- seq(0,R, length.out = 1000)
  # 
  # multiple_R <- integrate_theta(R_seq = R_seq,
  #                               alpha = alpha,
  #                               rconstraints = rconstraints,
  #                               Nupper = Nupper)
  # 
  # 
  # return( sum(multiple_R)/length(R_seq))
}



#check out how our integration changes with competition coeff
multiple_omega<-function(diagonal,off_diagonal, R, rconstraints=NULL, Nupper=NULL){
  
  alpha <- alpha<- diag(diagonal, ncol = 2, nrow = 2)
  inter_alphas <- seq(0, off_diagonal, 0.01)
  
  omegas <- sapply(inter_alphas, function(x, R, rconstraints, Nupper){
    
    alpha[1, 2] <- x
    alpha[2, 1] <- x
    
    saavedra_integration <- Omega_SA(alpha)
    analytical_intecration <- omega_analytical(alpha)
    numerical_integration <- integrate_radii(
      alpha = alpha,
      R =R,
      rconstraints = rconstraints,
      Nupper = Nupper )
    
    
    results <-  c( x,saavedra_integration, analytical_intecration, numerical_integration) %>% as.data.frame()
    
    return( results)
    
  },R=R, rconstraints =rconstraints, Nupper=Nupper   )
  
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

multiple_radii<- function(R, alpha, rconstraints=NULL, Nupper=NULL){
  R_seq<- seq(0.000001,R,length.out = 100) %>% as.list()
  
  values<- lapply(R_seq, function(x, alpha, rconstraints, Nupper){
    saavedra <- Omega_SA(alpha)
    theta_integration <- integrate_radii(alpha=alpha,
                                         R=x,
                                         rconstraints = rconstraints,
                                         Nupper= Nupper)
    
  results <- data.frame("R"=x, "Saavedras"=saavedra, "Theta"= theta_integration) 
  print(results)
  return(results)
    
    
  } ,alpha=alpha,
  rconstraints = rconstraints,
  Nupper = Nupper)
  
  calcs <- do.call(rbind, values) %>% as.data.frame()
  
  plot( calcs$R, calcs$Saavedras, pch=16, col="firebrick",
        ylab= "Omega",
        xlab= "Radius",
        ylim = c(0,1))
  points(calcs$R, calcs$Theta, pch=18, col="dodgerblue")
  
  return(tt)
  }

r_feasible<-function(alpha, rconstraints=NULL, Nupper=NULL, make_plot=FALSE){
  R_vals <- c(0.0001,1:100)
  
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
  
  return(r_sample)
}





rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(1, 1)
)
Nupper <- c(
  i = Inf,
  j = Inf
)







test<-multiple_omega( diagonal = 1,
                    off_diagonal = .9,
                    R = 1, 
                    rconstraints = rconstraints, 
                    Nupper =  Nupper)



test_radii <- multiple_radii(R = 100,alpha = alpha,rconstraints = rconstraints,Nupper = Nupper)






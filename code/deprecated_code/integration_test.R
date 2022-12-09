# integration tests

source("code/integration_toolbox.R")



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


#How the integation changes as R changes
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


#A visualization of the constraints 
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




#we make explicit the constraints of the growth rate of the two species
rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf)
)

#As well as their maximum abundances
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






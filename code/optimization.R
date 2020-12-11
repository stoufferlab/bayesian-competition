library(profvis)


profvis({
feasibility_theta <- function(theta_seq,alpha, R,rconstraints=NULL,Nupper=NULL){
  #we do not use solve to save time because posteriors are looooong
  inverse_alpha <- inverse_matrix(alpha)
  #we integrate over values of theta
  results <- sapply(theta_seq, function(theta, inverse_alpha, R, rconstraints, Nupper){
    #print(theta)
    ri <- R * cos(theta)
    rj <- R * sin(theta)
    r <- c(ri, rj)
    
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
      } else{
        #solve fo abundances
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
               N_good <- (N <= Nupper) 
               N_good <- all(N_good)
             } else{
               N_good <- TRUE
             }
          return(N_good)
        }
        
   }
  },inverse_alpha=inverse_alpha, R=R, rconstraints = rconstraints, Nupper=Nupper)
  
  
  return(results)
  
  
}

integrate_theta <-function( R_seq,alpha,rconstraints=NULL,Nupper=NULL ) {
  #we break down a circle into very small parts
  thetas <- seq(0 , 2*pi, length.out = 2000)
  
  results <- matrix(sapply(R_seq, function(R, alpha, rconstraints, Nupper){
    
    area<- feasibility_theta(theta_seq = thetas,
                             alpha = alpha,
                             R = R,
                             rconstraints = rconstraints,
                             Nupper= Nupper)
    
    # print(sum(area)/length(thetas))
    return(sum(area)/length(thetas))
  }, alpha=alpha,
  rconstraints = rconstraints,
  Nupper = Nupper))
  
  # results<- as.matrix(results)
  
  return(results)
  
}



rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(1, 1)
)


Nupper <- c(
  i = 1000,
  j = 1000
)

alpha <- matrix(c(0.026687122, 0.02897693 ,0.012431050, 0.03833583),nrow=2)

#ptm <- proc.time()
test1<- integrate_theta(alpha = alpha,
                        rconstraints=rconstraints,
                        Nupper= Nupper,
                        R_seq = seq(0,2000,1))

#proc.time() - ptm


})

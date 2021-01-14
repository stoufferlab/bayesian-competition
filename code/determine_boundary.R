require("sp")
source("code/integration_toolbox.R")

# R_max is the maximum radius calculated for a particular combo of constraints and alpha matrix
feasibility_boundary <- function(theta,
                                 alpha, 
                                 R_max,
                                 rconstraints=NULL,
                                 Nupper=NULL){
  # for every theta value we evaluate if ri and rj (determined by theta and R) along an R sequence
  N <- 500
  R_seq <- seq(0.0001, R_max, length.out = N) %>% as.list()
  R_boundary <- lapply(R_seq,
                       function(R,
                                theta,
                                alpha,
                                rconstraints=NULL,
                                Nupper=NULL){
                         
                         ri <- R * cos(theta)
                         rj <- R * sin(theta)
                         r <- c(ri, rj)
                         
                         #is this particular magnitude feasible
                         feasible <- check_feasibility(r = r,
                                                       alpha = alpha, 
                                                       rconstraints = rconstraints,
                                                       Nupper = Nupper)
                         
                         results <- data.frame( "theta" = theta, 
                                                "R_bound"= R,
                                                "ri"=ri,
                                                "rj"=rj ,
                                                "feasible"= feasible)
                         return(results)
                         
                       } ,theta = theta,
                       alpha = alpha,
                       rconstraints = rconstraints,
                       Nupper = Nupper)
  
                     R_bounds <- do.call(rbind, R_boundary)
                     
                    # #which values are unfeasible? 
                    #  unfeasible <- filter(R_bounds, feasible == 0)
                    #   
                    #  #if none, then return the maximum value of R
                    #  if(nrow(unfeasible) == 0){
                    #    Rb <- R_bounds[which(R_bounds$R_bound == max(R_bounds$R_bound) ), ]
                    #    bounds <- data.frame("theta"= theta, 
                    #                         "R_bound"= Rb[,"R_bound"],
                    #                         "ri"= Rb[,"ri"],
                    #                         "rj"= Rb[,"rj"])
                    #    return(bounds)
                    #  }else{
                    #    #if not return the minimum value of R where the growth rates are no longer feasible 
                    #    Rb <- unfeasible[which(unfeasible$R_bound == min(unfeasible$R_bound) ), ]
                    #    bounds <- data.frame("theta"= theta, 
                    #                         "R_bound"= Rb[,"R_bound"],
                    #                         "ri"= Rb[,"ri"],
                    #                         "rj"= Rb[,"rj"])
                    #    return(bounds)
                    #  }
                     
}

alpha_boundary <- function(alpha,
                           R_max,
                           rconstraints=NULL,
                           Nupper=NULL){
  # we add the shape delimited by the alphas alone
  # The column vectors of the alpha matrix define the boundary of the feasibility domain 
  # minus a tiny fraction because exactly in the boundary the feasibility is 0!
  m1 <- (alpha[2,2] / alpha[1,2]) - 1e-5 
  m2 <- alpha[2,1] / alpha[1,1]   - 1e-3
  
  values <- seq(0, R_max, length.out = 1000) %>% as.list()
  
  upper_limit <- lapply(values, function(ri,
                                             alpha,
                                             R_max, 
                                             rconstraints= rconstraints,
                                             Nupper= Nupper){
    rj <- ri*m1
    r<- c(ri, rj)
     feas <- check_feasibility(r = c(ri, rj),
                               alpha = alpha, 
                               rconstraints = rconstraints,
                               Nupper =  Nupper )
     if(feas){
       upper_limit<- data.frame("ri"=ri, "rj"=rj)
       return(upper_limit)
     }else{
       NULL
     }
    
    
  }, alpha=alpha,
  R_max =  R_max,
  rconstraints =  rconstraints,
  Nupper= Nupper)
  #now for the lower  defined by m2
  lower_limit <- lapply(values, function(ri,
                                         alpha,
                                         R_max, 
                                         rconstraints= rconstraints,
                                         Nupper= Nupper){
    rj <- ri*m2
    r<- c(ri, rj)
    feas <- check_feasibility(r = r,
                              alpha = alpha, 
                              rconstraints = rconstraints,
                              Nupper =  Nupper )
    if(feas){
      upper_limit<- data.frame("ri"=ri, "rj"=rj)
      return(upper_limit)
    }else{
      NULL
    }
    
    
  }, alpha=alpha,
  R_max =  R_max,
  rconstraints =  rconstraints,
  Nupper= Nupper)
  
  upper_points <- do.call(rbind, upper_limit)
  lower_points <- do.call(rbind, lower_limit)
}

feasibility_shape<-function( alpha, 
                             R_max,
                             rconstraints=NULL,
                             Nupper=NULL){
  N <- 500
  thetas <- seq(0 , 2*pi, length.out = N) %>% as.list()
  
  bound <- lapply(thetas, function(t,
                                    alpha,
                                    R_max, 
                                    rconstraints= rconstraints,
                                    Nupper= Nupper){
     bounded_points  <- feasibility_boundary(theta = t,
                                             alpha = alpha,
                                             R_max =  R_max,
                                             rconstraints =  rconstraints,
                                             Nupper= Nupper)
     return(bounded_points)
    
    
                                            
  }, alpha=alpha,
  R_max =  R_max,
  rconstraints =  rconstraints,
  Nupper= Nupper)
  
  shape <- do.call(rbind, bound)
  
  # we add the shape delimited by the alphas alone
  # The column vectors of the alpha matrix define the boundary of the feasibility domain 
  s1 <- alpha[2,2] / alpha[1,2]
  s2 <- alpha[2,1] / alpha[1,1]
  
  limit_rj <- shape[1,]
  ri_s1 <- seq(0, limit_rj[,"ri"], length.out = N/2)
  rj_s1 <- ri_s1 * s1
  
  upper_limit <- data.frame("ri"= ri_s1, "rj"=rj_s1)
  
  limit_ri <- filter(shape, ri == min(rj))
  ri_s2 <- seq(0, limit_rj[,"ri"], length.out = N/2)
  rj_s2 <- ri_s1 * s2
  
  lower_limit <- data.frame("ri"= ri_s2, "rj"=rj_s2)
  # All together, these coordinates describe the conical hull 
  
  coordinates <- shape[,c("ri","rj")] %>% rbind( upper_limit,lower_limit)
  
  return(coordinates)
}


alpha <- diag(2)
alpha[1,2] <- 0.3
alpha[2,1] <- 0.5
R_max <- 1000
rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf)
)

#And each species its maximum expected abundances
Nupper <- c(i = 200,
            j = 300)

r_feasible(alpha = alpha, rconstraints = rconstraints, Nupper = Nupper, R_max = R_max, make_plot = TRUE)




gg <- feasibility_shape(alpha = alpha , R_max = R_max , rconstraints = rconstraints, Nupper = Nupper)

points(gg$ri, gg$rj, pch=20, col= "dodgerblue")




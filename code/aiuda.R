
library(brms)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#We source everything known to human kind...


source("code/read_models.R")
source("code/feasibility_toolbox.R")
source("code/saavedra_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")
source("code/model_combo.R")



vero_model <- vero_lv_multispecies_poisson.rds
trcy_model <- trcy_bh_multispecies_poisson.rds


#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324

Ni_max<- 1e4
Nj_max<- 1e4

env <-  FALSE
bounded <- TRUE

alpha <- matrix( c(0.01006203,0.003556041, 0.07510000, 0.070489992), nrow = 2, ncol = 2, byrow = TRUE)

get_integration <- function(N_max, 
                            alpha,
                            bounded){
  
  
  #isaacs mean estimates
  gi<- 0.9641188
  si<- 0.9654804
  gj<- 0.4743974
  sj<- 0.9693324
  
  Ni_max<- 1e4
  Nj_max<- 1e4
  
  env <-  FALSE

  if(!bounded) {
    # Each model has its own constraints
    rconstraints <- list(
      lower = c(-Inf, -Inf),
      upper = c(Inf, Inf)
    )
    
  }else{
    # Each model has its own constraints
    rconstraints <- list(
      lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
      upper = c(vero_model$constraints[2], trcy_model$constraints[2])
    )  
  }
  
  #And each species its maximum expected abundances
  Nupper <- c(i = N_max,
              j = N_max)
  
  #Which determine the Radius for each alpha matrix used
  R_mean <- determine_radius(alpha = alpha,
                             Ni_max = N_max,
                             Nj_max = N_max)
  omega_SA_mean <- Omega_SA(alpha = alpha_mean)
  
  plot(0,0,
       xlim=c(-R_mean,R_mean),
       ylim=c(R_mean, R_max),
       type='n',
       xlab=expression(italic(r[i])),
       ylab=expression(italic(r[j]))
  )
  abline(h=0,lty='dashed',lwd=1.5)
  abline(v=0,lty='dashed',lwd=1.5)
  
  
  #Now we determine the feasibility shape with our Monte Carlo Integration
  integration_mean<- integrate_area(R_max = R_mean,
                                    alpha = alpha,
                                    rconstraints = rconstraints,
                                    Nupper = Nupper,
                                    n_samples = 1e4
  )
  
  
  #Now we determine the feasibility shape with our Monte Carlo Integration
  Omega_mean <- integration_mean$proportion
  #And also the coordinates of all the points that are feasible
  shape_mean <- integration_mean$coords
  # which tell us the bounds of the feasibility domain
  shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
  bounds_mean <- shape_bounds_mean$bounds
  # and also the area the area of the feasibility domain
  area_mean <- shape_bounds_mean$area
  
  
  #we store the values of coexistence using the point estimates
  mean_parameters_results <- data.frame(
    "Omega_saavedraa_mean"= omega_SA_mean,
    "Omega_mean"= Omega_mean,
    "area_mean"= area_mean,
    "R_mean"=R_mean,
    "N_max"= N_max)
  
  return(mean_parameters_results)
  
}

t1<- get_integration(N_max = 100, alpha = alpha,bounded = TRUE)


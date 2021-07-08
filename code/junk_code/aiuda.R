
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

check_point <- function(r,R_max,inv_alpha,rconstraints=NULL,Nupper=NULL){
  #returns NA if point is outside boundary
  #FALSE if it is inside boundary but infeasible
  #TRUE if it is inside boundary and feasible
  
  if(!check_radius_boundaries(r = r,
                              R_max = R_max)){
    #print("out of r boundaries")
    
    return(NA)
  }
  
  
  
  if(!check_r_boundaries(r = r,
                         rconstraints = rconstraints)){
    #  print("out of growth boundaries")
    points(r[1],r[2], pch=20, col = rethinking::col.alpha("#1f00c6",alpha=0.1))
    return(NA)
  }
  
  #solve fo abundances
  N  <- calculate_abundances(r = r,
                             inv_alpha = inv_alpha)
  
  if(!check_N_boundaries(N = N,
                         Nupper = Nupper)){
    # print("out of abundance boundaries")
    points(r[1],r[2], pch=20, col = rethinking::col.alpha("#c60044",alpha=0.1))
    return(NA)
  }
  
  N_feasible <- (N > 0)
  N_feasible <- all(N_feasible)
  
  return(N_feasible)
}

make_domain<- function(R,
                       alpha,
                       rconstraints,
                       Nupper){
  plot(0,0,
       xlim=c(-R, R),
       ylim=c(-R, R),
       type='n',
       xlab=expression(italic(r[i])),
       ylab=expression(italic(r[j]))
  )
  abline(h=0,lty='dashed',lwd=1.5)
  abline(v=0,lty='dashed',lwd=1.5)
  
  
  integration_mean<- integrate_area(R_max = R,alpha = alpha,rconstraints = rconstraints,Nupper = Nupper,desired_feasible = 1000,max_samples = 1e5)
  
  Omega_mean <- integration_mean$proportion
  
  
  #And also the coordinates of all the points that are feasible
  shape_mean <- integration_mean$coords
  #and the points that are unfeasible
  unfeasible_mean <- integration_mean$unfeasible
  # which tell us the bounds of the feasibility domain
  shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
  bounds_mean <- shape_bounds_mean$bounds
  # and also the area the area of the feasibility domain
  area_mean <- shape_bounds_mean$area
  
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=0.1)
  col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
  points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
  points(unfeasible_mean$ri, unfeasible_mean$rj, pch=20, col=col2)
  lines(bounds_mean$ri, bounds_mean$rj, col= "black", lwd=2)
  
  return(Omega_mean)
  
}


vero_model <- vero_rc_multispecies_poisson.rds
trcy_model <- trcy_rc_multispecies_poisson.rds


#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324

Ni_max<- 1e3
Nj_max<- 1e3

env <-  FALSE
bounded <- TRUE


alpha <- get_fixed_alphas(
  vero_model = vero_model,
  trcy_model = trcy_model,
  gi = gi,
  gj = gj,
  env = env)



rconstraints <- list(
  lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
  upper = c(vero_model$constraints[2], trcy_model$constraints[2])
)  

Nupper <- c(i = Ni_max,
            j = Nj_max)
#Which determine the Radius for each alpha matrix used
R  <- determine_radius(alpha = alpha,
                           Ni_max = Ni_max,
                           Nj_max = Nj_max)

make_domain(R = R,
            alpha = alpha,
            rconstraints = rconstraints,
            Nupper= Nupper)

# New proporitons ---------------------------------------------------------




get_boundary_r <-function(intraspecific_competition,
                           N_max,
                           lower,
                           upper){
  max_growth_rate <- intraspecific_competition * N_max
  
  if( intraspecific_competition > 0){

    bounds <- min(max_growth_rate, upper)
    return(bounds)
    
  }else{
    bounds <- max(max_growth_rate, lower)
   return(bounds)
  }
  
  return(bounds)
  
}


ri_bound <- get_boundary_r(intraspecific_competition = alpha[1,1],
                           N_max = Ni_max,
                           lower = rconstraints$lower[1],
                           upper = rconstraints$upper[1])


rj_bound <- get_boundary_r(intraspecific_competition = alpha[2,2],
                           N_max = Nj_max,
                           lower = rconstraints$lower[2],
                           upper = rconstraints$upper[2])





plot(0,0,
     xlim=c(-R,R),
     ylim=c(-R,R),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)


rect(xleft = 0,
     xright = ri_bound,
    ytop = rj_bound,
    ybottom = 0,
    col = rethinking::col.alpha("grey50", 0.5))



check_point <- function(r,R_max,inv_alpha,rconstraints=NULL,Nupper=NULL){
  #returns NA if point is outside boundary
  #FALSE if it is inside boundary but infeasible
  #TRUE if it is inside boundary and feasible
  
  if(!check_radius_boundaries(r = r,
                              R_max = R_max)){
    #print("out of r boundaries")
    
    return(NA)
  }
  
  
  
  if(!check_r_boundaries(r = r,
                         rconstraints = rconstraints)){
    #  print("out of growth boundaries")
    #points(r[1],r[2], pch=20, col = rethinking::col.alpha("#1f00c6",alpha=0.1))
    return(NA)
  }
  
  #solve fo abundances
  N  <- calculate_abundances(r = r,
                             inv_alpha = inv_alpha)
  
  if(!check_N_boundaries(N = N,
                         Nupper = Nupper)){
    # print("out of abundance boundaries")
    #points(r[1],r[2], pch=20, col = rethinking::col.alpha("#c60044",alpha=0.1))
    return(NA)
  }
  
  N_feasible <- (N > 0)
  N_feasible <- all(N_feasible)
  
  return(N_feasible)
}

integration_mean<- integrate_area(R_max = R,
                                  alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper,
                                  desired_feasible = 2000,
                                  max_samples = 3e5)


#And also the coordinates of all the points that are feasible
shape_mean <- integration_mean$coords

# which tell us the bounds of the feasibility domain
shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
bounds_mean <- shape_bounds_mean$bounds
# and also the area the area of the feasibility domain
area_mean <- shape_bounds_mean$area


col1 <- rethinking::col.alpha("mediumseagreen", alpha=0.1)
col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
lines(bounds_mean$ri, bounds_mean$rj, col= "black", lwd=2)


area_bounds <- ri_bound * rj_bound

proportion_area <- area_mean/ area_bounds

proportion_area

Omega_SA(alpha)



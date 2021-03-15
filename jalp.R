library(tidyverse)
source("code/model_toolbox.R")
source("code/determine_radius.R")
source("code/saavedra_toolbox.R")
source("code/read_models.R")
source("code/feasibility_toolbox.R")

rr <- results_sunny_bounded %>% filter(vero_model =="Lotka-Volterra") %>%
  filter(trcy_model =="Lotka-Volterra")

vero_model <- vero_lv_multispecies_poisson.rds
trcy_model <- trcy_lv_multispecies_poisson.rds

Ni_max <- 1e3
Nj_max <- 1e3

#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324

env <- FALSE
bounded <- TRUE



#we determine the alpha matrix for the mean parameter values
alpha_mean <- get_fixed_alphas(
  vero_model = vero_model,
  trcy_model = trcy_model,
  gi = gi,
  gj = gj,
  env = env)
#as well as r1 (vero's growth rate)
vero_growth_mean <- get_fixed_growth(
  model = vero_model,
  s = si,
  g = gi,
  env = env)
#and r2 (trcy's growth rate)
trcy_growth_mean <- get_fixed_growth(
  model = trcy_model,
  s = sj,
  g = gj,
  env = env)
#store them in a vector
r_mean <- c(vero_growth_mean,
            trcy_growth_mean)
print(r_mean)
#determine if there are boundaries on the growth rates
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
Nupper <- c(i = Ni_max,
            j = Nj_max)
#Which determine the Radius for each alpha matrix used
R_mean <- determine_radius(alpha = alpha_mean,
                           Ni_max = Ni_max,
                           Nj_max = Nj_max)
print(R_mean)
#Saaveda et al. estimation, to compare
omega_SA_mean <- Omega_SA(alpha = alpha_mean)

centroid_SA_mean <- r_centroid(alpha_mean)

theta_SA_mean <- theta(r_c = centroid_SA_mean,
                       r = r_mean)

feasibility_SA_mean <- test_feasibility_saavedra(alpha = alpha_mean,
                                                 r = r_mean)
#check if our growth rates are feasible, accordig to our constraints
feasiblity_mean <- check_point(r =r_mean,
                               R_max = R_mean,
                               inv_alpha =  inverse_matrix(alpha_mean),
                               rconstraints = rconstraints,
                               Nupper = Nupper)
feas_na <- is.na(feasiblity_mean)
#in this case if it is NA it means it is unfeasible
feasiblity_mean <-  ifelse(feas_na,0, feasiblity_mean)

#Now we determine the feasibility shape with our Monte Carlo Sampling
integration_mean<- integrate_area_converge(R_max = R_mean,
                                           alpha = alpha_mean,
                                           rconstraints = rconstraints,
                                           Nupper = Nupper,
                                           desired_feasible = 500,
                                           max_samples = 1e6
)
#which spits out the propotion of the area that is feasible, or the feasibility domain
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
#with the bounds we can then get the distance from the limit of our growth rates
distances_mean <- distance_from_limit(r=r_mean,
                                      shape = bounds_mean,
                                      feasibility = feasiblity_mean)



plot(0,0,
     xlim=c(-R_mean, R_mean),
     ylim=c(-R_mean, R_mean),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)


col1 <- rethinking::col.alpha("grey50", alpha=0.5)
points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
lines(bounds_mean$ri, bounds_mean$rj, col= "mediumseagreen", lwd=2)




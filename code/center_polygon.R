source("code/integration_toolbox.R")
require(polylabelr)


Nupper <- c(i = 100,
            j = 100)

rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = Inf, Inf)

alpha<- diag(2)
alpha[2,1]<- 0.3
alpha[1,2]<- 0.3

R_max <- 200

r_feasible(alpha = alpha,
           rconstraints = rconstraints,
           Nupper = Nupper,
           R_max = R_max,
           make_plot = TRUE)

shape <- feasibility_shape(alpha = alpha,
                          R_max = R_max,
                          rconstraints = rconstraints, 
                          Nupper = Nupper)

center <- poi(x = shape$ri,
    y = shape$rj)




source("code/integration_toolbox.R")
source("code/determine_boundary.R")
source("code/determine_radius.R")
require("polylabelr")



rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf))

alpha<- diag(2)
alpha[2,1]<- 0.3
alpha[1,2]<- 0.3


r_feasible(alpha = alpha,
           rconstraints = rconstraints,
           Nupper = NULL,
           R_max = 2,
           make_plot = TRUE)

shape <- feasibility_shape(alpha = alpha,
                          R_max = 2,
                          rconstraints = rconstraints, 
                          Nupper = NULL)



center <- poi(x = shape$ri,
    y = shape$rj)

center_distance <- shortest_distance(r= c(center$x, center$y),
                  shape = shape)


points(center$x, 
       center$y,
       pch=20 )

r<- c(1.5, 0.8)

growth_distance <- shortest_distance(r = r,
                                     shape = shape)


points(r[1], 
       r[2],
       pch=18 )

normalized_distance <- center_distance - growth_distance

normalized_distance


##########################################################

Ni_max <- 3
Nj_max <- 3

rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(1, Inf))


r<- c(1.1, 0.8)

Nupper <- c(i = Ni_max,
            j = Nj_max)

R_max<- determine_radius(alpha = alpha,
                         Ni_max = Ni_max,
                         Nj_max = Nj_max)


plot(0,0,
     xlim=c(-6,6),
     ylim=c(-6,6),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)


shape <- feasibility_shape(alpha = alpha,
                           R_max = R_max,
                           rconstraints = rconstraints, 
                           Nupper = Nupper)

lines(shape$ri, shape$rj)
points(r[1], r[2], pch=20)

#the center of the polygon
center <- poi(x = shape$ri,
              y = shape$rj)

points(center$x, center$y, pch=17)
#we get the shortest distance from the center of the polygon to an edge
center_distance <- shortest_distance(r= c(center$x, center$y),
                                     shape = shape)
#we get the shortest distance from our growth rates to an edge
growth_distance <- shortest_distance(r = r,
                                     shape = shape)


fixed_distance <- distance_from_limit(alpha = alpha,
                                      R_max =  R_max,
                                      rconstraints = rconstraints,
                                      Nupper =  Nupper,
                                      r =r )

results <- data.frame("center_distance" =center_distance,
                      "growth_distance"= growth_distance)

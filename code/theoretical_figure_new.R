source("code/feasibility_toolbox.R")

R<-2

rconstraints <- list(
  lower = c(-Inf, -1),
  upper = c(1, Inf))

alpha<- diag(2)
alpha[2,1]<- 0.3
alpha[1,2]<- 0.3

Ni_max <- 1e3
Nj_max <- 1e3

Nupper <- c(i = Ni_max,
            j = Nj_max)


r<- c(0.7,0.7)



make_projection<- function(R,
                           alpha,
                           rconstraints,
                           Nupper,
                           r){
  



integration_mean<- integrate_area(R_max = R,
                                           alpha = alpha,
                                           rconstraints = rconstraints,
                                           Nupper = Nupper,
                                           desired_feasible = 5000,
                                           max_samples = 1e6
)

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
distances_mean <- distance_from_limit(r=r,
                                      shape = bounds_mean,
                                      feasibility = 1)



plot(0,0,
     xlim=c(-R, R),
     ylim=c(-R, R),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)

col1 <- rethinking::col.alpha("grey50", alpha=0.1)
points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
lines(bounds_mean$ri, bounds_mean$rj, col= "mediumseagreen", lwd=2)
points(r[1],r[2], col="#e3004e", pch=18)

return(list(area_mean, distances_mean, integration_mean$proportion))

}


make_projection(R = R,alpha = alpha,rconstraints = rconstraints,Nupper = Nupper,r = r)
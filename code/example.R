Ni_max <- Inf
Nj_max <- Inf

rconstraints <- list(
  lower = c(-Inf, -1),
  upper = c(1, Inf))


Nupper <- c(i = Ni_max,
            j = Nj_max)
R_max<- 10


alpha <- matrix(c(0.011455927,0.09678921,  0.0004742236, 0.07965776),
                nrow = 2,
                ncol = 2,
                byrow = FALSE)

#alpha <- diag(2)
#alpha[1,2] <- 0.1

r <- c(0.7885806 , 35.40820)


plot(0,0,
     xlim=c(-R_max,R_max),
     ylim=c(-R_max,R_max),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)

norm_area<- integrate_area(R_max = R_max,
                           alpha = alpha,
                           rconstraints = rconstraints,
                           Nupper = Nupper,
                           n_samples = 1e5
)

proportion <- norm_area$proportion
shape <- norm_area$coords

shape_bounds <-determine_boundary_shape(shape)

area <- shape_bounds$area
bounds <- shape_bounds$bounds

lines(bounds$ri, bounds$rj)
points(r[1], r[2])

feasible <- check_point(r =r, R_max = R_max, inv_alpha = inverse_matrix(alpha),rconstraints = rconstraints,Nupper = Nupper)

distances<-distance_from_limit(r= r,
                    shape = bounds,
                    feasibility = feasible)

#the proportion of the feasibility domain
proportion
#the area
area
#are growth rates feasible
feasible
#distance form limit
distances$growth_distance



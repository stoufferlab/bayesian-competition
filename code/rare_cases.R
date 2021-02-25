

test_posterior <- test_2[order(decreasing = TRUE,test_2$distance_growth),]

rare <- test_2[2,]

 vero_model <- vero_lv_multispecies_poisson.rds
 trcy_model <- trcy_bh_multispecies_poisson.rds

 rconstraints <- list(
   lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
   upper = c(vero_model$constraints[2], trcy_model$constraints[2])
 )
 
 Nupper <- c(i = 1e4,
             j = 1e4)
 
 
alpha <- matrix(c(0.003247007,0.02957824, 0.00287961,0.04696119 ), nrow = 2, ncol = 2)
r <-c(0.6829131, 7.877205)
R <- 767.8425


r_feasible(alpha = alpha,
           rconstraints = rconstraints,
           Nupper = Nupper,
           R_max = R,
           make_plot = TRUE)


shape <- feasibility_shape(alpha = alpha,
                           R_max = R,
                           rconstraints = rconstraints, 
                           Nupper = Nupper)


distance <- shortest_distance(r = r,
                              shape = shape,
                              col = "dodgerblue")

plot(0,0,
     xlim=c(0,2),
     ylim=c(0,20),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)


lines(shape$ri, shape$rj)
points(r[1], r[2], pch=20)



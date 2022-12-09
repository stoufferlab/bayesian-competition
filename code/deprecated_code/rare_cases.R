

 vero_model <- vero_lv_multispecies_poisson.rds
 trcy_model <- trcy_bh_multispecies_poisson.rds

 rconstraints <- list(
   lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
   upper = c(vero_model$constraints[2], trcy_model$constraints[2])
 )
 
 Nupper <- c(i = 1e4,
             j = 1e4)
 
 
alpha <- matrix(c(0.003404525 ,0.02933614  ,0.0050805000 ,0.04467113), nrow = 2, ncol = 2)
r <-c(0.6746920 ,8.488632 )
R <- 744




shape <- determine_feasibility_shape(alpha = alpha,
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
points(r_mean[1], r_mean[2], pch=20)



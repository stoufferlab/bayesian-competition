

lv_bh <- mod %>% filter(vero_model == "Lotka-Volterra") %>% filter(trcy_model == "Beverton-Holt")

lv_bh <- lv_bh %>% filter( feasibility == TRUE)

ggplot(lv_bh) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = distance_growth,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  )


tt<- lv_bh[
  order( lv_bh$distance_growth, decreasing = TRUE ),
  ]
# +
#   geom_point(mapping = aes(x = Omega_mean,
#                            y = distance_growth_mean),
#              col = "goldenrod3") +
#   theme_alba +
#   scale_color_manual(values = c(col2, col1)) +
#   geom_abline(
#     intercept = 0,
#     slope = 0,
#     linetype = "dashed",
#     col = "grey50"
#   )


alpha <- matrix(c(0.011455927,0.09678921,  0.0004742236, 0.07965776),
                nrow = 2,
                ncol = 2,
                byrow = FALSE)

r <- c(0.7885806 , 35.40820)

R_max <- 1768.498

Ni_max <- 1e4
Nj_max <- 1e4

rconstraints <- list(
  lower = c(-Inf, -1),
  upper = c(1, Inf))


Nupper <- c(i = Ni_max,
            j = Nj_max)

R_max<- determine_radius(alpha = alpha,
                         Ni_max = Ni_max,
                         Nj_max = Nj_max)

plot(0,0,
     xlim=c(-1,1),
     ylim=c(0,40),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j]))
)
abline(h=0,lty='dashed',lwd=1.5)
abline(v=0,lty='dashed',lwd=1.5)

shape <- determine_feasibility_shape(alpha = alpha,
                           R_max = R_max,
                           rconstraints = rconstraints, 
                           Nupper = Nupper)

lines(shape$ri, shape$rj)
points(r[1], r[2], pch=20)



check_feasibility(r = r,
                  alpha = alpha,rconstraints = rconstraints,Nupper = Nupper)




r_feasible<-function(alpha, rconstraints=NULL, Nupper=NULL,R_max ,make_plot=FALSE){
  R_vals <- seq(0, R_max, length.out = 5000)
  
  #We sample values of R that are feasible to calculate their median, or the area in the center
  r_sample <- t(sapply(
    seq_len(5000),
    function(x,R_vals,alpha,rconstraints,Nupper){
      while(TRUE){
        R <- sample(R_vals,1)
        theta <- runif(1,0,2*pi)
        feas <- feasibility_theta(theta = theta,
                                  R= R,
                                  alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper)
        if(feas){
          return(c(R=R,theta=theta,ri=R*cos(theta),rj=R*sin(theta)))
        }
      }
    },
    R_vals=R_vals,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  )) %>% as.data.frame()
  
  
  #medianR <- Gmedian(r_sample[,c("ri","rj")])
  
  if(make_plot){
    plot(0,0,
         xlim=c(-R_max,R_max),
         ylim=c(-R_max,R_max),
         type='n',
         xlab=expression(italic(r[i])),
         ylab=expression(italic(r[j]))
    )
    abline(h=0,lty='dashed',lwd=1.5)
    abline(v=0,lty='dashed',lwd=1.5)
    apply(r_sample, MARGIN=1,
          function(x) {
            segments(0,0,x["ri"],x["rj"],col=grey(0.75))
          }
    )
  }
  # return(medianR)
  # return(r_sample)
}


r_feasible(alpha = alpha,rconstraints = rconstraints,Nupper = Nupper,R_max = 50,make_plot = TRUE)
points(r[1], r[2], pch=20)
shape <- determine_feasibility_shape(alpha = alpha,
                                     R_max = 50,
                                     rconstraints = rconstraints, 
                                     Nupper = Nupper)
lines(shape$ri, shape$rj)

##############################################################################

Omega_SA(alpha)
omega_analytical(alpha)


alpha<- diag(2)

r_feasible(alpha = alpha,rconstraints = rconstraints,Nupper = Nupper,R_max = 2,make_plot = TRUE)



R_seq <- seq(0.1,2,0.1)
area <-c()


for(i in 1:length(R_seq)){
  tt <- feasibility_theta(theta_seq = thetas, 
                          alpha = alpha, 
                          R= R_seq[i], 
                          rconstraints = rconstraints, 
                          Nupper = Nupper)
  area <- c(area, sum(tt)/length(tt))

  
  
}



multiple_R <- hcubature(f=integrate_theta,
                        lowerLimit  = c(0,0), 
                        upperLimit  = c(2,2),
                        vectorInterface = TRUE,
                        tol = 1e-3,
                        alpha = alpha,
                        rconstraints= rconstraints,
                        Nupper= Nupper)














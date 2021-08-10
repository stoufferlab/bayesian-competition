source("code/feasibility_toolbox.R")
source("code/figure_label.R")
source("code/determine_radius.R")

area_species_alone <- function(alpha,
                               Nupper,
                               rconstraints){
  
  
  ri_bound <- get_boundary_r(intraspecific_competition = alpha[1,1],
                             N_max = Nupper[1],
                             lower = rconstraints$lower[1],
                             upper = rconstraints$upper[1])
  
  
  rj_bound <- get_boundary_r(intraspecific_competition = alpha[2,2],
                             N_max = Nupper[2],
                             lower = rconstraints$lower[2],
                             upper = rconstraints$upper[2])
  
  area_bounds <-abs(  ri_bound * rj_bound)
  
  return(c(ri_bound,
           rj_bound))
  
  
}

structural_stability_wrapper_figure <- function(
                                         alpha,
                                         rconstraints,
                                         Ni_max,
                                         Nj_max,
                                         r,
                                         label){
  
  plot(0,0,
       xlim=c(-3, 3),
       ylim=c(-3, 3),
       type='n',
       xlab="",
       ylab="",
       cex.lab=1.5
  )
  title(main = label, adj=0, line = 0.5, font.main=1, cex.main = 1.5)
  abline(h=0,lty='dashed',lwd=1.5)
  abline(v=0,lty='dashed',lwd=1.5)
  
  points(r[1], r[2], pch=19,)
  
  
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  
  R <- determine_radius(alpha = alpha,
                        Ni_max = N_max,
                        Nj_max = N_max)
  
  
  #first we check if our observed growth rates are feasible
  feasiblity <- check_point(r =r,
                            R_max = R,
                            inv_alpha =  inverse_matrix(alpha),
                            rconstraints = rconstraints,
                            Nupper = Nupper)
  feas_na <- is.na(feasiblity)
  #in this case if it is NA it means it is unfeasible
  feasiblity <-  ifelse(feas_na,0, feasiblity)
  
  #we do a mcm integration of the feasible area
  integration<- integrate_area(R_max = R,
                               alpha = alpha,
                               rconstraints = rconstraints,
                               Nupper = Nupper,
                               desired_feasible = 4000,
                               max_samples = 5e5
  )
  
  area_alone <- area_species_alone(alpha = alpha,
                                   Nupper = Nupper,
                                   rconstraints = rconstraints)


  
  aa <- integration$coords
  
  
  shape_bounds <- determine_boundary_shape(shape = integration$coords)
  bounds <- shape_bounds$bounds
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=1)
 
  
  points(aa$ri, aa$rj, col= col1, pch=19)
  
  
  points(r[1], r[2], pch=19)
  
  
  
  #rect(xleft = 0,
   #    xright = area_alone[1],
    #   ytop = area_alone[2],
     #  ybottom = 0,
      # border = NA,
       #col = rethinking::col.alpha("grey50", alpha = 0.3))
  # and also the area the area of the feasibility domain
  area <- shape_bounds$area
  #with the bounds we can then get the distance from the limit of our growth rates
  distances <- distance_from_limit(r=r,
                                   shape = bounds,
                                   feasibility = feasiblity)
  distance_from_edge <- distances$growth_distance
  
  #but also we get the proportion of things inside the convex hull
  convex_mean <- calculate_convex(shape = bounds,
                                  unfeasible = integration$unfeasible)
  
  
  proportion <- area/ (area_alone[1] * area_alone[2])
  
  
  results <- data.frame("area_feasible"= area,
                        "area_alone"= area_alone[1] * area_alone[2],
                        "proportion"= proportion,
                        "convex"= convex_mean,
                        "feasibility"= feasiblity,
                        "distance"= distance_from_edge,
                        "detection"= distances$detection)
  
  return(results)
  
  
  
}


pdf(file="../bayesian_competition_ms/rho_delta.pdf",width = 8, height = 6)
layout.matrix <- matrix(c(1,2,3,
                          4,4,4),
                        nrow = 2, ncol = 3, byrow = T)
layout(mat = layout.matrix, heights=c(1,1,1), widths = c(1,1))



alpha<- diag(2)
alpha[1,1] <- 1
alpha[2,2] <- 1
alpha[1,2]<- 0.5
alpha[2,1]<- 0.5


rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))

one <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    r = c(1, 1),
    label = "A)",
    Ni_max = 1.5,
    Nj_max = 1.5
  )


rho_one <- one$area_feasible / (2.5 * 2.5)
delta_one <- one$distance

rect(xleft = 0,
    xright = 2.5,
   ytop = 2.5,
    ybottom = 0,
 border = NA,
col = rethinking::col.alpha("grey50", alpha = 0.3))

# two ---------------------------------------------------------------------


alpha<- diag(2)
alpha[1,1] <- 1
alpha[2,2] <- 1
alpha[1,2]<- 0.7
alpha[2,1]<- 0.9


rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))

two <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    r = c(2, 0.7),
    label = "B)",
    Ni_max = 1.5,
    Nj_max = 1.5
  )


rho_two <- two$area_feasible / (2.5 * 2.5)
delta_two <- two$distance

rect(xleft = 0,
     xright = 2.5,
     ytop = 2.5,
     ybottom = 0,
     border = NA,
     col = rethinking::col.alpha("grey50", alpha = 0.3))




# three -------------------------------------------------------------------



alpha<- diag(2)
alpha[1,1] <- 1
alpha[2,2] <- 1
alpha[1,2]<- 0.2
alpha[2,1]<- 0.3


rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))

three <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    r = c(2, 0.7),
    label = "C)",
    Ni_max = 1.7,
    Nj_max = 2
  )


rho_three <- three$area_feasible / (1 * 2.3)
delta_three <- three$distance

rect(xleft = 0,
     xright = 1,
     ytop = 2.3,
     ybottom = 0,
     border = NA,
     col = rethinking::col.alpha("grey50", alpha = 0.3))




mtext(expression("Growth rate"~italic(r[i])), side = 1, outer = TRUE, line = -24,cex = 1)


mtext(expression("Growth rate"~italic(r[j])~"                     "), side = 2, outer = TRUE, line = -2.5, cex = 1, adj = 1)



# final -------------------------------------------------------------------

par(mar=c(4,18,2,18))

plot(0,0,
     xlim=c(0, 1.5),
     ylim=c(-1, 1),
     type='n',
     xlab="",
     ylab="",
     cex.lab=1.5
)
title(main = "D)", adj=0, line = 0.5, font.main=1, cex.main = 1.5)


abline(0,0, lwd =1.2, col= "grey50")

text( x= rho_one, y= delta_one, "A", cex = 1.5)
text( x= rho_two, y= delta_two, "B", cex = 1.5)
text( x= rho_three, y= delta_three, "C", cex = 1.5)


mtext(expression("Relative coexistence ratio,"~rho), side = 1, outer = TRUE, line = -1,cex = 1)

mtext(expression("               Distance from the edge,"~delta), side = 2, outer = TRUE, line = -16,cex = 1, adj = 0)


dev.off()

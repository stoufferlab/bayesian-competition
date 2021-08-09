source("code/feasibility_toolbox.R")
source("code/figure_label.R")
source("code/determine_radius.R")



structural_stability_wrapper_figure <- function(
                                         alpha,
                                         rconstraints,
                                         N_max,
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
  
  
  Ni_max <- N_max
  Nj_max <-N_max
  
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
  
  aa <- integration$coords
  
  
  shape_bounds <- determine_boundary_shape(shape = integration$coords)
  bounds <- shape_bounds$bounds
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=1)
 
  
  points(aa$ri, aa$rj, col= col1, pch=19)
  
  
  points(r[1], r[2], pch=19)
  
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
  
  
  area_alone <- area_species_alone(alpha = alpha,
                                   Nupper = Nupper,
                                   rconstraints = rconstraints)
  
  proportion <- area/ area_alone
  
  
  results <- data.frame("area_feasible"= area,
                        "area_alone"= area_alone,
                        "proportion"= proportion,
                        "convex"= convex_mean,
                        "feasibility"= feasiblity,
                        "distance"= distance_from_edge,
                        "detection"= distances$detection)
  
  return(results)
  
  
  
}




alpha<- diag(2)
alpha[1,2]<- 0.5
alpha[2,1]<- 0.5


rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))

one <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    N_max = 1.5,
    r = c(1, 1),
    label = "A)"
  )

rect(xleft = 0,
     xright = 2.5,
     ytop = 2.5,
     ybottom = 0,
     border = NA,
     col = rethinking::col.alpha("grey50", alpha = 0.3))

rho_one <-    one/ (2.5*2.5)


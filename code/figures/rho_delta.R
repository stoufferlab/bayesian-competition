
require(here)

# defines some functions used below
source(here("code/lib/feasibility_utils.R"))

structural_stability_wrapper_figure <- function(
                                         alpha,
                                         rconstraints,
                                         giNi_max,
                                         gjNj_max,
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
  
  # points(r[1], r[2], pch=19,)
  
  
  gN_max <- c(i = giNi_max,
              j = gjNj_max)
  
  R <- maximum_radius(alpha = alpha,
                      gN_max = gN_max)
  
  #we do a mcm integration of the feasible area
  integration<- integrate_area(
    alpha = alpha,
    R_max = R,
    rconstraints = rconstraints,
    gN_max = gN_max,
    desired_feasible = 5000,
    max_samples = 50000
  )
  
  area_alone <- area_species_alone(alpha = alpha,
                                   gN_max = gN_max,
                                   rconstraints = rconstraints)


  
  aa <- integration$coords
  
  
  shape_bounds <- determine_boundary_shape(shape = integration$coords)
  bounds <- shape_bounds$bounds
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=1)
 
  # DEBUG: replace with convex hull
  # bcfd <- apply(feasible_pts, 1, within_N_boundaries, inv_alpha=inv_alpha, gN_max=gN_max)
  # bcfd_pts <- feasible_pts[bcfd,]
  bcfd_hull <- grDevices::chull(aa)
  bcfd_hull <- c(bcfd_hull, bcfd_hull[1])

  polygon(
    aa[bcfd_hull,],
    border = NA,
    lty=0,
    col = scales::alpha("mediumseagreen", alpha = 1)
  )
  # points(aa$ri, aa$rj, col= col1, pch=19)
  
  # observed growth rates
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
                                   shape = bounds)
  distance_from_edge <- distances$growth_distance
  if(!is_feasible(r, inverse_matrix(alpha)))
    distance_from_edge <- -1 * distance_from_edge
  
  # #but also we get the proportion of things inside the convex hull
  # convex_mean <- calculate_convex(shape = bounds,
  #                                 unfeasible = integration$unfeasible)
  
  
  proportion <- area/ (area_alone[1] * area_alone[2])
  
  
  results <- data.frame("area_feasible"= area,
                        "area_alone"= area_alone[1] * area_alone[2],
                        "proportion"= proportion,
                        # "convex"= convex_mean,
                        "distance"= distance_from_edge,
                        "detection"= distances$detection)
  
  return(results)
  
  
  
}


setEPS()
postscript(file=here("figures/rho_delta.ps"),width = 8, height = 6)
layout.matrix <- matrix(c(1,2,3,
                          4,4,4),
                        nrow = 2, ncol = 3, byrow = T)
layout(mat = layout.matrix, heights=c(1,1,1), widths = c(1,1))

rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf)
)

alpha<- diag(2)
alpha[1,1] <- 1
alpha[2,2] <- 1
alpha[1,2]<- 0.5
alpha[2,1]<- 0.7

rvals <- c(1, 1.2)
gN_max <- c(1.5,1.5)



one <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    r = rvals,
    label = "A)",
    giNi_max = gN_max[1],
    gjNj_max = gN_max[2]
  )


rho_one <- one$area_feasible / (alpha[1,1] * gN_max[1] + alpha[2,2] * gN_max[2])
delta_one <- one$distance

rect(
  xleft = 0,
  xright = alpha[1,1] * gN_max[1],
  ytop = alpha[2,2] * gN_max[2],
  ybottom = 0,
  border = NA,
  col = rethinking::col.alpha("grey50", alpha = 0.3)
)

# two ---------------------------------------------------------------------


alpha<- diag(2)
alpha[1,1] <- 1
alpha[2,2] <- 1
alpha[1,2]<- 0.7
alpha[2,1]<- 0.9

gN_max <- c(1.5,1.5)

# rconstraints <- list(
#   lower = c(-1, -1),
#   upper = c(Inf, Inf))

two <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    r = c(2, 0.7),
    label = "B)",
    giNi_max = gN_max[1],
    gjNj_max = gN_max[2]
  )


rho_two <- two$area_feasible / (alpha[1,1] * gN_max[1] + alpha[2,2] * gN_max[2])
delta_two <- two$distance

rect(
  xleft = 0,
  xright = alpha[1,1] * gN_max[1],
  ytop = alpha[2,2] * gN_max[2],
  ybottom = 0,
  border = NA,
  col = rethinking::col.alpha("grey50", alpha = 0.3)
)



# three -------------------------------------------------------------------



alpha<- diag(2)
alpha[1,1] <- 1
alpha[2,2] <- 1
alpha[1,2]<- 0.4
alpha[2,1]<- 0.2


# rconstraints <- list(
#   lower = c(-1, -1),
#   upper = c(Inf, Inf))

rvals <- c(1.5, 1.3)

gN_max <- c(2.5,2.5)

three <-
  structural_stability_wrapper_figure(
    alpha = alpha,
    rconstraints = rconstraints,
    r = rvals,
    label = "C)",
    giNi_max = gN_max[1],
    gjNj_max = gN_max[2]
  )


rho_three <- three$area_feasible / (alpha[1,1] * gN_max[1] + alpha[2,2] * gN_max[2])
delta_three <- three$distance

rect(
  xleft = 0,
  xright = alpha[1,1] * gN_max[1],
  ytop = alpha[2,2] * gN_max[2],
  ybottom = 0,
  border = NA,
  col = rethinking::col.alpha("grey50", alpha = 0.3)
)

mtext(expression("Vital rate,"~italic(r[i])), side = 1, outer = TRUE, line = -25,cex = 1)


mtext(expression("Vital rate,"~italic(r[j])~"                     "), side = 2, outer = TRUE, line = -2.1, cex = 1, adj = 1)



# final -------------------------------------------------------------------

par(mar=c(4,18,2,18))

plot(NA,NA,
     xlim=c(0.2, 5),
     ylim=c(-1, 1),
     type='n',
     xlab="",
     ylab="",
     log='x',
     cex.lab=1.5
)
title(main = "D)", adj=0, line = 0.5, font.main=1, cex.main = 1.5)


abline(h=0, lwd =1.5) #, col= "grey50")

text( x= rho_one, y= delta_one, "A", cex = 1.5)
text( x= rho_two, y= delta_two, "B", cex = 1.5)
text( x= rho_three, y= delta_three, "C", cex = 1.5)


mtext(expression("Relative coexistence ratio,"~rho), side = 1, outer = TRUE, line = -1,cex = 1)

mtext(expression("               Distance from the edge,"~delta), side = 2, outer = TRUE, line = -16,cex = 1, adj = 0)


dev.off()

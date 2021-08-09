source("code/feasibility_toolbox.R")
source("code/figure_label.R")
source("code/determine_radius.R")

pdf(file="../bayesian_competition_ms/feasibility_domain.pdf",width = 10, height = 3)


#xlab=expression(italic(r[i]))
#ylab=expression(italic(r[j]))

make_figure_area <- function(alpha,
                            N_max,
                             rconstraints,
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
  
  
  
  Ni_max <- N_max
  Nj_max <-N_max
  
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  
   R <- determine_radius(alpha = alpha,
                         Ni_max = N_max,
                         Nj_max = N_max)
  
  
  integration_mean<- integrate_area(R_max = R,
                                    alpha = alpha,
                                    rconstraints = rconstraints,
                                    Nupper = Nupper,
                                    desired_feasible = 40000,
                                    max_samples = 9e5)
  
  #And also the coordinates of all the points that are feasible
  shape_mean <- integration_mean$coords
  
  # which tell us the bounds of the feasibility domain
  shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
  bounds_mean <- shape_bounds_mean$bounds
  # and also the area the area of the feasibility domain
  area_mean <- shape_bounds_mean$area
  
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=1)
  col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
  points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
  
  
  #fig_label(text = label, region = "figure", pos  ="top", cex=1.5)
}


layout.matrix <- matrix(c(1,2,3,4),
                        nrow = 1, ncol = 4, byrow = T)
layout(mat = layout.matrix, heights=c(1,1,1), widths = c(1,1))

# feasibility domain ------------------------------------------------------


alpha<- diag(2)
alpha[1,2]<- 0.5
alpha[2,1]<- 0.5



rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf))

make_figure_area(alpha = alpha,
                 N_max = 6,
                 rconstraints = rconstraints,
                 label = "A) Feasibility domain")





# abundance ---------------------------------------------------------------

make_figure_area <- function(alpha,
                             N_max,
                             rconstraints,
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
  
  
  
  Ni_max <- N_max
  Nj_max <-N_max
  
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  
  R <- determine_radius(alpha = alpha,
                        Ni_max = N_max,
                        Nj_max = N_max)
  
  
  integration_mean<- integrate_area(R_max = R,
                                    alpha = alpha,
                                    rconstraints = rconstraints,
                                    Nupper = Nupper,
                                    desired_feasible = 8000,
                                    max_samples = 5e5)
  
  #And also the coordinates of all the points that are feasible
  shape_mean <- integration_mean$coords
  
  # which tell us the bounds of the feasibility domain
  shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
  bounds_mean <- shape_bounds_mean$bounds
  # and also the area the area of the feasibility domain
  area_mean <- shape_bounds_mean$area
  
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=1)
  col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
  #points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
  
  
  #fig_label(text = label, region = "figure", pos  ="top", cex=1.5)
}
rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf))

check_point <- function(r,R_max,inv_alpha,rconstraints=NULL,Nupper=NULL){
  #returns NA if point is outside boundary
  #FALSE if it is inside boundary but infeasible
  #TRUE if it is inside boundary and feasible
  
  if(!check_radius_boundaries(r = r,
                              R_max = R_max)){
  
    
    if(r[1]>0 && r[2]>0){
      points(r[1],r[2], pch=20, col = rethinking::col.alpha("#c60044",alpha=1))
    } 
     #points(r[1],r[2], pch=20, col = rethinking::col.alpha("#c60044",alpha=1))
    
    return(NA)
  }
  
  
  
  if(!check_r_boundaries(r = r,
                         rconstraints = rconstraints)){
    #  print("out of growth boundaries")
     # points(r[1],r[2], pch=20, col = rethinking::col.alpha("#1f00c6",alpha=0.1))
    return(NA)
  }
  
  #solve fo abundances
  N  <- calculate_abundances(r = r,
                             inv_alpha = inv_alpha)
  
  if(!check_N_boundaries(N = N,
                         Nupper = Nupper)){
    # print("out of abundance boundaries")
    if(r[1]>0 && r[2]>0){
      points(r[1],r[2], pch=20, col = rethinking::col.alpha("#c60044",alpha=1))
    } 
    #p
    
    
    return(NA)
  }
  
  N_feasible <- (N > 0)
  N_feasible <- all(N_feasible)
  
  return(N_feasible)
}



make_figure_area(alpha = alpha,
                 N_max = 1.5,
                 rconstraints = rconstraints,
                 label = "B) Abundance constraints")

# model -------------------------------------------------------------------

#par(  mar=c(4,4,2,1))

check_point <- function(r,R_max,inv_alpha,rconstraints=NULL,Nupper=NULL){
  #returns NA if point is outside boundary
  #FALSE if it is inside boundary but infeasible
  #TRUE if it is inside boundary and feasible
  
  if(!check_radius_boundaries(r = r,
                              R_max = R_max)){
    
    #  points(r[1],r[2], pch=20, col = rethinking::col.alpha("#1f00c6",alpha=1))
    
    
    
    return(NA)
  }
  
  
  
  if(!check_r_boundaries(r = r,
                         rconstraints = rconstraints)){
    #  print("out of growth boundaries")
    #points(r[1],r[2], pch=20, col = rethinking::col.alpha("#1f00c6",alpha=0.1))
    return(NA)
  }
  
  #solve fo abundances
  N  <- calculate_abundances(r = r,
                             inv_alpha = inv_alpha)
  
  if(!check_N_boundaries(N = N,
                         Nupper = Nupper)){
    # print("out of abundance boundaries")
    
    #  points(r[1],r[2], pch=20, col = rethinking::col.alpha("#c60044",alpha=1))
    return(NA)
  }
  
  N_feasible <- (N > 0)
  N_feasible <- all(N_feasible)
  
  return(N_feasible)
}

rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))
make_figure_area(alpha = alpha,
                 N_max = 1.5,
                 rconstraints = rconstraints,
                 label = "C) Model-based constraints")


rect(xleft = -4,
     xright = -1,
     ytop = 4,
     ybottom = -4,
     border = NA,
     col = rethinking::col.alpha("#1f00c6", alpha = 1))

rect(xleft = -4,
     xright = 4,
     ytop = -1,
     ybottom = -4,
     border = NA,
     col = rethinking::col.alpha("#1f00c6", alpha = 1))


# 
# rconstraints <- list(
#   lower = c(-1, -Inf),
#   upper = c(Inf, 1))
# make_figure_area(alpha = alpha,
#                  N_max = 1.5,
#                  rconstraints = rconstraints,
#                  label = "D")
# 
# rect(xleft = -4,
#      xright = 4,
#      ytop = 4,
#      ybottom = 1,
#      border = NA,
#      col = rethinking::col.alpha("#1f00c6", alpha = 1))
# 
# rect(xleft = -4,
#      xright = -1,
#      ytop = 4,
#      ybottom = -4,
#      border = NA,
#      col = rethinking::col.alpha("#1f00c6", alpha = 1))

# both -------------------------------------------------------------------


make_figure_area <- function(alpha,
                             N_max,
                             rconstraints,
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
  
  
  
  Ni_max <- N_max
  Nj_max <-N_max
  
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  
  R <- determine_radius(alpha = alpha,
                        Ni_max = N_max,
                        Nj_max = N_max)
  
  
  integration_mean<- integrate_area(R_max = R,
                                    alpha = alpha,
                                    rconstraints = rconstraints,
                                    Nupper = Nupper,
                                    desired_feasible = 8000,
                                    max_samples = 5e5)
  
  #And also the coordinates of all the points that are feasible
  shape_mean <- integration_mean$coords
  
  # which tell us the bounds of the feasibility domain
  shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
  bounds_mean <- shape_bounds_mean$bounds
  # and also the area the area of the feasibility domain
  area_mean <- shape_bounds_mean$area
  
  
  col1 <- rethinking::col.alpha("mediumseagreen", alpha=1)
  col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
  points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
  
  
  #fig_label(text = label, region = "figure", pos  ="top", cex=1.5)
}


check_point <- function(r,R_max,inv_alpha,rconstraints=NULL,Nupper=NULL){
  #returns NA if point is outside boundary
  #FALSE if it is inside boundary but infeasible
  #TRUE if it is inside boundary and feasible
  
  if(!check_radius_boundaries(r = r,
                              R_max = R_max)){
    
    if(r[1]>0 && r[2]>0){
    #  points(r[1],r[2], pch=20, col = rethinking::col.alpha("darkmagenta",alpha=1))
    } 
    
  
    return(NA)
  }
  
  
  
  if(!check_r_boundaries(r = r,
                         rconstraints = rconstraints)){
  #  points(r[1],r[2], pch=20, col = rethinking::col.alpha("darkmagenta",alpha=0.1))
    return(NA)
  }
  
  #solve fo abundances
  N  <- calculate_abundances(r = r,
                             inv_alpha = inv_alpha)
  
  if(!check_N_boundaries(N = N,
                         Nupper = Nupper)){
    if(r[1]>0 && r[2]>0){
     # points(r[1],r[2], pch=20, col = rethinking::col.alpha("darkmagenta",alpha=1))
    } 
    return(NA)
  }
  
  N_feasible <- (N > 0)
  N_feasible <- all(N_feasible)
  
  return(N_feasible)
}



rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))

make_figure_area(alpha = alpha,
                 N_max = 1.5,
                 rconstraints = rconstraints,
                 label = "D) BCFD")



mtext(expression("Growth rate"~italic(r[i])), side = 1, outer = TRUE, line = -1,cex = 1.3)


mtext(expression("Growth rate"~italic(r[j])), side = 2, outer = TRUE, line = -2.5, cex = 1.3)

dev.off()

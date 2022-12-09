
require(here)

# defines some functions used below
source(here("code/lib/feasibility_utils.R"))

pdf(file=here("figures/feasibility_domain.pdf"),width = 10, height = 3)

#xlab=expression(italic(r[i]))
#ylab=expression(italic(r[j]))

make_figure_area <- function(label){
  plot(0,0,
       xlim=c(-3, 3),
       ylim=c(-3, 3),
       type='n',
       xlab="",
       ylab="",
       cex.lab=1.5
  )
  # title(main = label, adj=0, line = 0.5, font.main=1, cex.main = 1.5)
  abline(h=0,lty='dashed',lwd=1.5)
  abline(v=0,lty='dashed',lwd=1.5)
  
  # #And also the coordinates of all the points that are feasible
  # shape_mean <- integration_mean$coords
  
  # # # which tell us the bounds of the feasibility domain
  # # shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
  # # bounds_mean <- shape_bounds_mean$bounds
  # # # and also the area the area of the feasibility domain
  # # area_mean <- shape_bounds_mean$area
  
  # col1 <- scales::alpha("mediumseagreen", alpha=1)
  # col2 <- scales::alpha("grey50", alpha = 0.1)
  # points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
  
  
  # fig_label(text = label, region = "figure", pos  ="top", cex=1.5)
}


layout.matrix <- matrix(c(1,2,3,4),
                        nrow = 1, ncol = 4, byrow = T)
layout(mat = layout.matrix, heights=c(1,1,1), widths = c(1,1))

# biologically constrained feasibility domain ------------------------------------------------------

alpha<- diag(2)
alpha[1,2]<- 0.25
alpha[2,1]<- 0.25

inv_alpha <- solve(alpha)

rconstraints <- list(
  lower = c(-Inf, -Inf),
  upper = c(Inf, Inf))

gN_max <- c(1, 2)

R_max <- maximum_radius(
  alpha = alpha,
  gN_max = c(10,10)
)

integration_mean<- integrate_area(
  alpha = alpha,
  R_max = R_max,
  rconstraints = rconstraints,
  gN_max = c(10,10),
  desired_feasible = 50000,
  max_samples = 100000
)

feasible_pts <- integration_mean$coords
excluded_pts <- integration_mean$unfeasible

# feasibility

make_figure_area()
mtext("Mathematical", side = 3, line=2)
mtext("feasibility domain", side = 3, line=0.8)
# mtext("Feasibility domain", side = 3, line=1.4)
mtext("                         - ", side = 3, line=1.2, cex=2)

apply(
  rbind(feasible_pts,excluded_pts),
  1,
  function(x){
    if(is_feasible(x, inv_alpha))
      points(
        x[1],
        x[2],
        pch=20,
        # col=scales::alpha("#c60044", alpha=1)
        col=scales::alpha("mediumseagreen", alpha=1)
      )
  }
)

# abundance

rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf))

make_figure_area()
mtext("Abundance constraints", side = 3, line=1.4)
mtext("                            -", side = 3, line=1.2, cex=2)

# polygon(
#   c(-4,-4,-1,-1,4,4,-4),
#   c(-4,4,4,-1,-1,-4,-4),
#   border = NA,
#   lty=0,
#   col = scales::alpha("#1f00c6", alpha = 1)
# )

apply(
  rbind(feasible_pts,excluded_pts),
  1,
  function(x){
    if(is_feasible(x, inv_alpha) && !within_N_boundaries(x, inv_alpha, gN_max))
      points(
        x[1],
        x[2],
        pch=20,
        col=scales::alpha("#c60044", alpha=1)
      )
  }
)

# rect(xleft = -1,
#      xright = 4,
#      ytop = -1,
#      ybottom = -4,
#      border = NA,
#      lty=0,
#      col = scales::alpha("#1f00c6", alpha = 0.1))


# # 
# # rconstraints <- list(
# #   lower = c(-1, -Inf),
# #   upper = c(Inf, 1))
# # make_figure_area(alpha = alpha,
# #                  N_max = 1.5,
# #                  rconstraints = rconstraints,
# #                  label = "D")
# # 
# # rect(xleft = -4,
# #      xright = 4,
# #      ytop = 4,
# #      ybottom = 1,
# #      border = NA,
# #      col = scales::alpha("#1f00c6", alpha = 1))
# # 
# # rect(xleft = -4,
# #      xright = -1,
# #      ytop = 4,
# #      ybottom = -4,
# #      border = NA,
# #      col = scales::alpha("#1f00c6", alpha = 1))

# model-based

make_figure_area()
mtext("Model-based constraints", side = 3, line=1.4)

apply(
  rbind(feasible_pts,excluded_pts),
  1,
  function(x){
    if(!within_r_boundaries(x, rconstraints))
      points(
        x[1],
        x[2],
        pch=20,
        col=scales::alpha("#1f00c6", alpha=1)
      )
  }
)

# bcfd in the house!

make_figure_area()

# mtext("A)", side = 3, at=c(-1.3,0), xpd=NA)
mtext("Biologically constrained", side = 3, line=2)
mtext("feasibility domain", side = 3, line=0.8)
mtext("=                            ", side = 3, line=1.2, cex=2)

# # if we want to plot the domain using points
# apply(
#   rbind(feasible_pts),
#   1,
#   function(x){
#     if(within_N_boundaries(x, inv_alpha, Nupper))
#       points(
#         x[1],
#         x[2],
#         pch=20,
#         # col=scales::alpha("#c60044", alpha=1)
#         col=scales::alpha("mediumseagreen", alpha=1)
#       )
#   }
# )

bcfd <- apply(feasible_pts, 1, within_N_boundaries, inv_alpha=inv_alpha, gN_max=gN_max)
bcfd_pts <- feasible_pts[bcfd,]
bcfd_hull <- grDevices::chull(bcfd_pts)
bcfd_hull <- c(bcfd_hull, bcfd_hull[1])

polygon(
  bcfd_pts[bcfd_hull,],
  border = NA,
  lty=0,
  col = scales::alpha("mediumseagreen", alpha = 1)
)

mtext(expression("Vital rate,"~italic(r[i])), side = 1, outer = TRUE, line = -2,cex = 1.1)
mtext(expression("Vital rate,"~italic(r[j])), side = 2, outer = TRUE, line = -2., cex = 1.1)

dev.off()

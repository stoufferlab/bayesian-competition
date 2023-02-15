
require(here)

# defines some functions used below
source(here("code/lib/feasibility_utils.R"))

pdf(file=here("figures/feasibility_domain.pdf"),width = 10, height = 3)

#xlab=expression(italic(r[i]))
#ylab=expression(italic(r[j]))

make_figure_area <- function(label){
  plot(NA,NA,
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
  mtext(expression("Vital rate,"~italic(r[i])), side = 1, line = 3,cex = 1.1)
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
  upper = c(Inf, Inf)
)

gN_max <- c(10, 10)

npts <- 5000

# sample feasible equilibria
N_feas <- sapply(gN_max, function(x,npts) runif(npts,0,x), npts=npts)

# determine the growth rates the correspond to these
r_feas <- N_feas %*% t(alpha)

# mathematical feasibility

make_figure_area()
mtext("Mathematical", side = 3, line=2)
mtext("feasibility domain", side = 3, line=0.7)
# mtext("Feasibility domain", side = 3, line=1.4)
mtext("                         - ", side = 3, line=1.2, cex=2)

# add in the feasible growth rates
points(
  r_feas[,1],
  r_feas[,2],
  pch=20,
  col=scales::alpha("mediumseagreen", alpha=1)
)

# abundance constraints

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

# impose new abundance constraints
gN_max <- c(1,2)

outside_abundance_constraints <- apply(
  N_feas,
  1,
  function(x,gN_max){any(x>gN_max)},
  gN_max=gN_max
)

# plot points outside the abundance constraints
points(
  r_feas[outside_abundance_constraints,1],
  r_feas[outside_abundance_constraints,2],
  pch=20,
  col=scales::alpha("#c60044", alpha=1)
)

# model-based constraints

rconstraints <- list(
  lower = c(-1, -1),
  upper = c(Inf, Inf)
)

make_figure_area()
mtext("Model-based constraints", side = 3, line=1.4)

r_null <- cbind(runif(npts,-5,5),runif(npts,-5,5))

r_bad <- r_null[apply(r_null, 1, function(x) x[1] <= -1 || x[2] <= -1),]

# plot points outside the vital rate constraints
points(
  r_bad[,1],
  r_bad[,2],
  pch=20,
  col=scales::alpha("#1f00c6", alpha=1)
)

# bcfd in the house!

make_figure_area()

mtext("Biologically constrained", side = 3, line=2)
mtext("feasibility domain", side = 3, line=0.7)
mtext("=                            ", side = 3, line=1.2, cex=2)

# sample feasible equilibria
r_feas <- integrate_area(
  alpha=alpha,
  rconstraints=rconstraints,
  gN_max=gN_max,
  npts
)$coords

bcfd_hull <- grDevices::chull(r_feas)
bcfd_hull <- c(bcfd_hull, bcfd_hull[1])
polygon(
  r_feas[bcfd_hull,],
  border = 'black',
  lty=0,
  col = scales::alpha("mediumseagreen", alpha = 1)
)

# bcfd <- apply(feasible_pts, 1, within_N_boundaries, inv_alpha=inv_alpha, gN_max=gN_max)
# bcfd_pts <- feasible_pts[bcfd,]
# bcfd_hull <- grDevices::chull(bcfd_pts)
# bcfd_hull <- c(bcfd_hull, bcfd_hull[1])

# polygon(
#   bcfd_pts[bcfd_hull,],
#   border = NA,
#   lty=0,
#   col = scales::alpha("mediumseagreen", alpha = 1)
# )

# mtext(expression("Vital rate,"~italic(r[i])), side = 1, outer = TRUE, line = -2,cex = 1.1)
mtext(expression("Vital rate,"~italic(r[j])), side = 2, outer = TRUE, line = -2., cex = 1.1)

dev.off()

source("code/feasibility_toolbox.R")
source("code/figure_label.R")
source("code/determine_radius.R")

pdf(file="../bayesian_competition_ms/posterior.pdf",width = 7, height = 7)



#par(oma=c(0,0,0,0),  mar=c(4,5,2,0.5))



layout.matrix <- matrix(c(1,2,3,
                          4,5,5), nrow = 2, ncol = 3, byrow = T)
layout(mat = layout.matrix, heights=c(1,1), widths = c(1,1,1)) 

par(pty="s", mar=c(4,4,2,1))

alpha<- diag(2)
alpha[1,2]<-0.5
alpha[2,1]<-0.5

r<- c(0.6, 0.6)

R<- 2.2

plot(0,0,
     xlim=c(0, R),
     ylim=c(0, R),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j])),
     cex.lab=1.5
)
#abline(h=0,lty='dashed',lwd=1.5)
#abline(v=0,lty='dashed',lwd=1.5)


rconstraints <- list(
  lower = c(-Inf, -1),
  upper = c(1, Inf))


Ni_max <-2
Nj_max <-2
Nupper <- c(i = Ni_max,
            j = Nj_max)


ri_bound <- get_boundary_r(intraspecific_competition = alpha[1,1],
                           N_max = Ni_max,
                           lower = rconstraints$lower[1],
                           upper = rconstraints$upper[1])


rj_bound <- get_boundary_r(intraspecific_competition = alpha[2,2],
                           N_max = Nj_max,
                           lower = rconstraints$lower[2],
                           upper = rconstraints$upper[2])


rect(xleft = 0,
     xright = ri_bound,
     ytop = rj_bound,
     ybottom = 0,
     border = NA,
     col = rethinking::col.alpha("grey50", 0.5))


R_det <- determine_radius(alpha = alpha,Ni_max = Ni_max,Nj_max = Nj_max
)

integration_mean<- integrate_area(R_max = R_det,
                                  alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper,
                                  desired_feasible = 5000,
                                  max_samples = 3e5)

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
points(r[1], r[2], pch=18, col="black", cex=2)

fig_label(text = "A", region = "figure", pos  ="topleft", cex=1.5)


p1<-structural_stability_wrapper(R= R_det,
                             alpha =  alpha,
                             rconstraints = rconstraints,
                             Nupper = Nupper,
                             r = r)

# B -----------------------------------------------------------------------

alpha<- diag(2)
alpha[1,2]<-0.7
alpha[2,2]<-0.7
alpha[2,1]<-0.9
alpha[2,2]<- 0.9
r<- c(0.6, .3)


plot(0,0,
     xlim=c(0, R),
     ylim=c(0, R),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j])),
     cex.lab=1.5
)
#abline(h=0,lty='dashed',lwd=1.5)
#abline(v=0,lty='dashed',lwd=1.5)


rconstraints <- list(
  lower = c(-Inf, -1),
  upper = c(1, Inf))


Ni_max <-2
Nj_max <-2
Nupper <- c(i = Ni_max,
            j = Nj_max)


ri_bound <- get_boundary_r(intraspecific_competition = alpha[1,1],
                           N_max = Ni_max,
                           lower = rconstraints$lower[1],
                           upper = rconstraints$upper[1])


rj_bound <- get_boundary_r(intraspecific_competition = alpha[2,2],
                           N_max = Nj_max,
                           lower = rconstraints$lower[2],
                           upper = rconstraints$upper[2])


rect(xleft = 0,
     xright = ri_bound,
     ytop = rj_bound,
     ybottom = 0,
     border = NA,
     col = rethinking::col.alpha("grey50", 0.5))


R_det <- determine_radius(alpha = alpha,Ni_max = Ni_max,Nj_max = Nj_max
)

integration_mean<- integrate_area(R_max = R_det,
                                  alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper,
                                  desired_feasible = 5000,
                                  max_samples = 3e5)


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
points(r[1], r[2], pch=17, col="black", cex=1.5)

fig_label(text = "B", region = "figure", pos  ="topleft", cex=1.5)


p2<-structural_stability_wrapper(R= R_det,
                             alpha =  alpha,
                             rconstraints = rconstraints,
                             Nupper = Nupper,
                             r = r)

# C -----------------------------------------------------------------------



alpha<- diag(2)
alpha[1,1]<- 0.7
alpha[2,2]<- 0.7
alpha[1,2]<-0.2
alpha[2,1]<-0.2

r<- c(0.9, 2.2)


plot(0,0,
     xlim=c(0, R),
     ylim=c(0, R),
     type='n',
     xlab=expression(italic(r[i])),
     ylab=expression(italic(r[j])),
     cex.lab=1.5
)
#abline(h=0,lty='dashed',lwd=1.5)
#abline(v=0,lty='dashed',lwd=1.5)


rconstraints <- list(
  lower = c(-Inf, -1),
  upper = c(1, Inf))


Ni_max <-2
Nj_max <-2
Nupper <- c(i = Ni_max,
            j = Nj_max)


ri_bound <- get_boundary_r(intraspecific_competition = alpha[1,1],
                           N_max = Ni_max,
                           lower = rconstraints$lower[1],
                           upper = rconstraints$upper[1])


rj_bound <- get_boundary_r(intraspecific_competition = alpha[2,2],
                           N_max = Nj_max,
                           lower = rconstraints$lower[2],
                           upper = rconstraints$upper[2])


rect(xleft = 0,
     xright = ri_bound,
     ytop = rj_bound,
     ybottom = 0,
     border = NA,
     col = rethinking::col.alpha("grey50", 0.5))


R_det <- determine_radius(alpha = alpha,Ni_max = Ni_max,Nj_max = Nj_max
)

integration_mean<- integrate_area(R_max = R_det,
                                  alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper,
                                  desired_feasible = 8000,
                                  max_samples = 3e5)

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
points(r[1], r[2], pch=15, col="black", cex=1.5)

fig_label(text = "C", region = "figure", pos  ="topleft", cex=1.5)


p3<-structural_stability_wrapper(R= R_det,
                             alpha =  alpha,
                             rconstraints = rconstraints,
                             Nupper = Nupper,
                             r = r)


# push --------------------------------------------------------------------
#par(pty="m", mar=c(4,4,6,1))

plot(0,0,
     xlim=c(1 ,2),
     ylim=c(0,2),
     type='n',
     ylab= "Multi-species feasible area",
     xlab= "Single-species feasible area",
     cex.lab=1.5
)


points(p1$area_alone, p1$area_feasible, pch=18, cex=2)
points(p2$area_alone, p2$area_feasible, pch=17, cex=1.5)
points(p3$area_alone, p3$area_feasible, pch=15, cex=1.5)


fig_label(text = "D", region = "figure", pos  ="topleft", cex=1.5)


# domain ------------------------------------------------------------------

par(pty="m", mar=c(4,8,2,1))

plot(0,0,
     xlim=c(0 ,1),
     ylim=c(-1,1),
     type='n',
     ylab= "Distance from the edge, D",
     xlab= "Relative coexistence area",
     cex.lab=1.5
)

abline(0, 0, lty=2, lwd=2, col="grey50")




points(p1$proportion, p1$distance, pch=18, cex=2)
points(p2$proportion, p2$distance, pch=17, cex=1.5)
points(p3$proportion, p3$distance, pch=15, cex=1.5)

fig_label(text="Coexistence", region = "plot", pos = "topleft", cex=1.5)
fig_label(text="Competitive exclusion", region = "plot", pos = "bottomleft", cex=1.5)


fig_label(text = "E", region = "figure", pos  ="topleft", cex=1.5)

dev.off()


source("code/contrained_toolbox.R")
#We create an arbitrary matrix
alpha <- diag(2)
alpha[1,2] <- -0.80
alpha[2,1] <- -0.223

# Lets assume there are no constraints on the growht rates
rconstraints <- list(
  lower = c( -Inf, -Inf),
  upper = c( Inf, Inf)
)

# are there upper bounds on N that we do not consider biologically realistic?
Nupper <- c(
  i = Inf,
  j = Inf
)
#We calculate and plot omega and its center, for an arbitrary R
p<-feasibility_wrapper( R = 100,
                     alpha = alpha,
                     rconstraints = rconstraints,
                     Nupper = Nupper,
                     plot=TRUE)
#Lets say this are our calculated growth rates
pp<-c(26,80)
points(pp[1],pp[2], pch=20, col="firebrick4")
#We check if they are feasible
check_feasibility(r = pp,
                  rconstraints = rconstraints,
                  Nupper = Nupper,
                  alpha = alpha )

#Now we change the constraints, lets say sp.1 is from a BH model and sp 2 from LV
rconstraints <- list(
  lower = c( -1, -Inf),
  upper = c( Inf, 1)
)


#and re do the analysis
feasibility_wrapper( R = 100,
                     alpha = alpha,
                     rconstraints = rconstraints,
                     Nupper = Nupper,
                     plot = TRUE)

pp<-c(26,80)
points(pp[1],pp[2], pch=20, col="firebrick4")
#They are no longer feasible!
check_feasibility(r = pp,
                  rconstraints = rconstraints,
                  Nupper = Nupper,
                  alpha = alpha )

########################### What about constraints on N#############################
rconstraints <- list(
  lower = c( -Inf, -Inf),
  upper = c( Inf, Inf)
)

# are there upper bounds on N that we do not consider biologically realistic?
Nupper <- c(
  i = Inf,
  j = Inf
)
#We calculate and plot omega and its center, for an arbitrary R
p<-feasibility_wrapper( R = 100,
                        alpha = alpha,
                        rconstraints = rconstraints,
                        Nupper = Nupper,
                        plot=TRUE)
#Lets say this are our calculated growth rates
pp<-c(26,80)
points(pp[1],pp[2], pch=20, col="firebrick4")
#We check if they are feasible
check_feasibility(r = pp,
                  rconstraints = rconstraints,
                  Nupper = Nupper,
                  alpha = alpha )





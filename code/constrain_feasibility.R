# load libraries ----------------------------------------------------------
library(tidyverse)
library(mvtnorm)
#This is the code Chuliant sent me regarding constrains and the feasibility domain
# functions ---------------------------------------------------------------
calculate_omega <- function(alpha) {
  S <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha)
  #not sure why he calculates these
  m <- matrix(0, S, 1)
  a <- matrix(0, S, 1)
  b <- matrix(Inf, S, 1)
  
  d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
  #this is the normalization he does so omega goes from 0 to 0.5
 d[1]^(1/S)
  #we used to return
  #return(d[1])
}

#This function samples inside the feasibility domain with the inequality constrains to see what percentage of the feasibility meets the constrains
calculate_constrained_proportion <- function(alpha, constraints){
  Nsample <- 1000
  # generate a uniform sample of growth rates per each constrain
  r <- constraints %>% 
    map(~runif(Nsample, min = .[1], max = .[2])) %>% 
    unlist() %>% 
    matrix(ncol = Nsample, byrow = TRUE) #every column is a sample (Nsample) and every row is per constrain
  # test whether the constrained r is feasible or not
  feasibility <- c()
  #for every sample we calculate the feasibility using the sampled growth rates and the interaction matrix
  for(i in 1:Nsample){
    N <- solve(alpha, -r[,i])
    feasibility[i] <- if_else(sum(N<0) == 0, 1, 0)
  }
  #chuliang returned the normalized feasibility, but lets just return the mean
  #mean(feasibility)^(1/nrow(alpha))
  return(mean(feasibility))
}


# illustration by Chuliant ------------------------------------------------------------

#generate a random interaction matrix
set.seed(1010)
alpha <- runif(3^2, min = -1, max = 0) %>% 
  matrix(ncol=3)
diag(alpha) <- -1

# generate a constraint on r
constraints <- list(c(-1,1), c(0,1), c(0,2)) 

# calculate the size of the full feasibility domain
calculate_omega(alpha)
# calculate the proportion of the constrained r inside feasibility domain
calculate_constrained_proportion(alpha, constraints)
# calculate the size of the constrained feasibility domain
calculate_omega(alpha) * calculate_constrained_proportion(alpha, constraints)


# Applying it to  our constrains ------------------------------------------------------------
# an actual matrix using mean paramter values for Beverton Holt both species
alpha <- matrix(c(0.026687122, 0.02897693 ,0.012431050, 0.03833583),nrow=2)

calculate_omega(alpha = alpha)
#what happens if they have the same constraints? shouldn't the proportion be 1?
x<- seq(10, 100000, 10)

shared_omega<-vector("numeric", length = length(x))
omegas<-calculate_omega(alpha)


for(i in 1:length(x)){
  constraints <- list(c(-x[i],x[i]), c(-x[i],x[i]))
  shared_omega[i]<- calculate_constrained_proportion(alpha = alpha,constraints = constraints)
  omega_prime[i] <- omegas * shared_omega[i]
  
}


#but it is not! it is 0! this does not make sense
calculate_constrained_proportion(alpha = alpha,constraints = constraints)

calculate_omega(alpha) * calculate_constrained_proportion(alpha, constraints)

# now we want to know how the proportion changes as x goes to infinty
#constraints <- list(c(-1,Inf), c(Inf,1))
# 
# x <- seq(0,100,1)
# proportion_que <- c()
# for(i in 1:length(x)){
#   constraints <- list(c(-1,x[i]), c(-x[i],1))
#   proportion  <- calculate_constrained_proportion(alpha = alpha,constraints = constraints)
#   proportion_que <- c(proportion_que, proportion)
# }
# 






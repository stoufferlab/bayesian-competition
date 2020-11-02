library(cubature)
library(Gmedian)


###Functions to calculate the feasibility domain, its center, and if our calculated growth rates are feasible############


#Checks if a particular combination of growth rates,given a Radius, is feasible with the interaction matrix and constraints
feasibility <- function(ri,alpha,R,branch=c("plus","minus"),rconstraints=NULL,Nupper=NULL){
  # growth rates are on a circle of radius R so given one we can determine the other
  rj <- switch(branch,
               plus  = sqrt(R^2 - ri^2),
               minus = -sqrt(R^2 - ri^2)
  )
  r <- c(ri,rj)
  
  # solve for the equilibrium given the interactions and the growth rate vector
  N <- solve(alpha) %*% r
  
  # check if N corresponds to feasibile equilibrium
  N_feasible <- (N > 0)
  
  # check whether N are within constraints
  if(!is.null(Nupper)){
    N_good <- (N <= Nupper)
  }else{
    N_good <- rep(TRUE,length(N))
  }
  
  # check whether growth rates are within constraints
  if(!is.null(rconstraints)){
    r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
  }else{
    r_good <- rep(TRUE,length(r))
  }
  
  return(prod(N_feasible, N_good, r_good))
}

#perform integration using the function above to determine the size of the feasibility domain
omega <- function(R,alpha,rconstraints=NULL,Nupper=NULL){
  
  
  
  # calculate the "plus" branch of the integral for the alpha matrix
  omega_plus <- adaptIntegrate(
    feasibility,
    lower=c(-R),
    upper=c(R),
    alpha=alpha,
    R=R,
    branch="plus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  # calculate the "minus" branch of the integral for the alpha matrix
  omega_minus <- adaptIntegrate(
    feasibility,
    lower=c(-R),
    upper=c(R),
    alpha=alpha,
    R=R,
    branch="minus",
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  # the total integral is the sum of the two
  omega_alpha <- omega_plus$integral + omega_minus$integral
  
  # because of the limits of the integral on r, we can normalize directly
  omega <- omega_alpha / (4*R)
  
  return(omega)
}

#Visualize how omega changes with different values of Radius (R_vals) given constraints
omega_radius<-function(R_vals, alpha,rconstraints=NULL, Nupper=NULL){

  
  # calculate the proportional feasibility domain across the different lengths of R
  omega_vals <- sapply(
    R_vals,
    omega,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  
  plot(R_vals, omega_vals, type='l', xlab=expression(italic(R)), ylab=expression(Omega))
  return(data.frame("R"=R_vals, "Omegas"= omega_vals, "Mean_omega"= mean(omega_vals)))
}

#Integrate the area of omega, R is the maximum radius expected
omega_integrate<-function(R, alpha, rconstraints=NULL, Nupper=NULL){
 
   omega_area<- adaptIntegrate(
    omega,
    lower=0,
    upper=R,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  )
  #we normalize it
  return(omega_area$integral/100)
}

# sample values of R that are feasible and get the center of the area, plot if you want
r_feasible<-function(R,alpha, rconstraints=NULL,Nupper=NULL, make_plot=FALSE){
  R_vals <- seq(0.001,R, 0.01)
  
  r_sample <- t(sapply(
    seq_len(100),
    function(x,R_vals,alpha,rconstraints,Nupper){
      while(TRUE){
        #you sample a value of R
        R <- sample(R_vals,1)

        theta <- runif(1,0,2*pi)

        ri <- R*cos(theta)
        branch <- ifelse(sin(theta)>0,"plus","minus")
        feas <- feasibility(ri,alpha,R,branch,rconstraints,Nupper)
        if(feas){
          return(c(R=R,theta=theta,ri=ri,rj=R*sin(theta)))
        }
      }
    },
    R_vals=R_vals,
    alpha=alpha,
    rconstraints=rconstraints,
    Nupper=Nupper
  ))
  
  if(make_plot){
    plot(0,0,
         xlim=c(-range(R_vals)[2],range(R_vals)[2]),
         ylim=c(-range(R_vals)[2],range(R_vals)[2]),
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
    # compute and plot the median vector (which appears to end up a bit short for some reason)
    medianR <- Gmedian(r_sample[,c("ri","rj")])
    points(medianR[1],medianR[2],pch=20,col="aquamarine4")
    
  }
  
 
  return(r_sample)
}


#Function to determine radius given maximum abundances N and an alpha matrix
determine_radius<-function(N, alpha){
  r<-alpha %*% N
  R<- sqrt(r[1]^2 + r[2]^2)
  return(R)
}
#Returns Omega and the center of the area
feasibility_wrapper<-function(R, alpha, rconstraints=NULL, Nupper=NULL, make_plot=FALSE){
  
  area <- omega_integrate (    R = R,
                               alpha = alpha,
                               rconstraints = rconstraints,
                               Nupper = Nupper)
  #We sample the values of R that are feasible
  growth_values<- r_feasible(R = R,
                             alpha = alpha,
                             rconstraints = rconstraints,
                             Nupper = Nupper,
                             make_plot = make_plot)
 
  #We get the median of the growth vectors to calculate the center of the area
  
  center<- Gmedian(growth_values[,c("ri","rj")])

  return(list("Omega"=area, "Center"= center))
  
}

#Check if our growth rates, r,  are feasible given alpha matrix and the constraints
check_feasibility<-function(r,alpha,rconstraints=NULL, Nupper=NULL){
  # solve for the equilibrium given the interactions and the growth rate vector
  N <- solve(alpha) %*% r
  
  # check if N corresponds to feasibile equilibrium
  N_feasible <- (N > 0)
  
  # check whether N are within constraints
  if(!is.null(Nupper)){
    N_good <- (N <= Nupper)
  }else{
    N_good <- rep(TRUE,length(N))
  }
  
  # check whether growth rates are within constraints
  if(!is.null(rconstraints)){
    r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
  }else{
    r_good <- rep(TRUE,length(r))
  }
  
  return(prod(N_feasible, N_good, r_good))
  
}

#Calculate the distance from the center 
calculate_distance<-function(center, r){
  dist<- sqrt((r[1]- center[1])^2 + (r[2]- center[2])^2)
  return(dist)
}



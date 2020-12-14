require(Rsolnp)


#to maximize f(x) is equivalent to minimizing -f(x)
# so this is a function that gives you -R^2
radius<- function(N, alpha, Ni_max, Nj_max){
  #Where N[1] is Ni and N[2] is Nj
  r1 = ( (N[1]*alpha[1,1]) + (N[2]*alpha[1,2]))^2
  r2 = ( (N[1]*alpha[1,2]) + (N[2]*alpha[2,2]))^2
  # remember that this actually gives you -R^2
  R = - r1 - r2
  return(R)
}

#Note that the main and constraint functions must take the exact same arguments, irrespective of whether they are used by all of them.
inequalities <- function(N,alpha, Ni_max, Nj_max) { 
  z1=  N[1] - Ni_max
  z2=  N[2] - Nj_max
  return(c(z1, z2))
}

#wrapper function that takes maximum values of abundances and the alpha matrix
#returns the value of R optimized
determine_radius<- function(alpha, Ni_max, Nj_max){
  # Ni/Nj has to be equal or smaller than Ni_max
  upper_constraints <- c(0,0)
  # Ni/Nj has to be equal or greater than 0
  lower_constraints <- c(-Ni_max,-Nj_max)
  #initial params at the center not at the boundary
  N0 = c(Ni_max/2,Nj_max/2)
  
  
  solution <- solnp(fun = radius,
                    pars = N0,
                    ineqfun = inequalities,
                    ineqUB = upper_constraints,
                    ineqLB = lower_constraints,
                    alpha = alpha,
                    Ni_max = Ni_max,
                    Nj_max = Nj_max)
  
  convergence <- solution$convergence
  
  if (convergence == 0) {
    #we extract R
    n <- length(solution$values)
    value  <- solution$values[2]
    R <- sqrt(-value)
    return(R)
  }else{
    warning("No convergence")
  }
  
  
}




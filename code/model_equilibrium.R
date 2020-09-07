
#monoculture equilibriums for different models
bh_equilibrium <-function(s,g,lambda, alphaii) {


     a <- (lambda*g)
     b <- 1 - ( (1-g)*s)
    equilibrium <- ((a/b)  -1) / (alphaii*g)

    return( equilibrium)
}


rc_equilibrium <-function(s,g,lambda, alphaii) {
  
  
  a <- (lambda*g)
  b <- 1 - ( (1-g)*s)
  equilibrium <- ( log(b/a) ) / -(alphaii*g)

  return( equilibrium)
}




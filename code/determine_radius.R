#determining the radius

determine_radius <-function(alpha, Ni_max, Nj_max){
  
  R1 <- sqrt( (alpha[1,1]*Ni_max)^2 + (alpha[2,1]*Ni_max)^2 )
  R2 <- sqrt( (alpha[1,2]*Nj_max)^2 + (alpha[2,2]*Nj_max)^2 )
  R3 <- sqrt ( (alpha[1,1]*Ni_max + alpha[1,2]*Nj_max)^2 + (alpha[2,1]*Ni_max +           alpha[2,2]*Nj_max)^2  )
  
  Ni <- (- Nj_max) * ( ((alpha[1,1]*alpha[1,2]) + (alpha[2,1]*alpha[2,2])) /(alpha[1,1]^2 + alpha[2,1]^2 )  )
  Nj <- (- Ni_max) * ( ((alpha[1,1]*alpha[1,2]) + (alpha[2,1]*alpha[2,2])) /(alpha[1,2]^2 + alpha[2,2]^2 )  )
  
  R4 <- sqrt ( (alpha[1,1]*Ni + alpha[1,2]*Nj_max)^2 + (alpha[2,1]*Ni + alpha[2,2]*Nj_max)^2  )
  R5 <- sqrt ( (alpha[1,1]*Ni_max + alpha[1,2]*Nj)^2 + (alpha[2,1]*Ni_max + alpha[2,2]*Nj)^2  )
  
  values <- c( R1, R2, R3, R4, R5)
  #print(values)
  return(max(values))
}
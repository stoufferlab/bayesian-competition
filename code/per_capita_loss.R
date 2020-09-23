#per capita loss functions of different models

lv_loss<-function(alpha, N){
  proportion_remaining <- (1-alpha*N)
  per_capita <- (1-proportion_remaining)/N
  
   return(per_capita)
}

rc_loss<-function(alpha, N){
  proportion_remaining<- (exp(-alpha*N))
  
  per_capita <- (1 -proportion_remaining) /N
 
  return(per_capita)
}

bh_loss<-function(alpha, N){
  proportion_remaining <- 1/(1+alpha*N)
  
  per_capita <- (1 -proportion_remaining) /N
  
  return(per_capita)
}

hs_loss<-function(alpha, N, beta){
  proportion_remaining <- 1/(1+alpha*N)^beta
  
  per_capita <- (1 -proportion_remaining) /N
  
  return(per_capita)
}



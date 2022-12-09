

inverse_matrix<-function(alpha){
   
   deteminant_alpha <- (alpha[1,1] * alpha[2,2]) - (alpha[2,1] * alpha[1,2])
   
 inverse_det  <- 1 / deteminant_alpha
  
  
  adjugate <- matrix( NA,nrow=2,ncol=2)                  
  adjugate[1,1] <- alpha[2,2]
  adjugate[2,2] <- alpha[1,1]
  
  adjugate[1,2] <- -alpha[1,2]
  adjugate[2,1] <- -alpha[2,1]
  
  
  ii <- inverse_det * adjugate
  
  return(ii)
  
  
}
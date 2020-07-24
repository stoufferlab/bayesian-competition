#get posterior distribution and add a column of growth rates:



#function to generate a posterior growth rate (and an environmental growth rate) for each point of the posterior, for one model. exp_param is binary that tells it if to take into consideration if the model has an exponential form (1=yes, 2=no)
posterior_parameters<-function(model, fun, s ,g, exp_param){
 post        <-posterior_samples(model)
 lambda      <- post$b_lambda_Intercept
 lambda_env <-  lambda + post$b_lambda_env

 if( isTRUE(exp_param)){
   b           <- post$b_b_Intercept
   b_env       <- post$b_b_env_Intercept
   growth     <- fun(s,g,lambda,b)
   env_growth <- fun(s,g,lambda,b_env)
 }else{
   growth     <- fun(s,g,lambda)
   env_growth <- fun(s,g,lambda_env)
 }
 
 all_posterior <- cbind(post,growth,env_growth)
 return(all_posterior)
}

#function to generate an alpha matrix based on one row of the posterior for each species. env is a binary that tells it if to take into consideration the environmental variables
alpha_matrix <- function(vero_row, trcy_row, gi,gj, env){
  # vero_row <- as.list(vero_row)
  # trcy_row <- as.list(trcy_row)
  
  if(env){
  alpha11 <- vero_row$b_alphaii_Intercept + vero_row$b_alphaii_env
  alpha21 <- trcy_row$b_alphaij_Intercept + trcy_row$b_alphaij_env
  alpha12 <- vero_row$b_alphaij_Intercept + vero_row$b_alphaij_env
  alpha22 <- trcy_row$b_alphaii_Intercept + trcy_row$b_alphaii_env
  alpha   <-matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2)
  alpha   <-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")

}else{
  alpha11 <- vero_row$b_alphaii_Intercept
  alpha21 <- trcy_row$b_alphaij_Intercept
  alpha12 <- vero_row$b_alphaij_Intercept
  alpha22 <- trcy_row$b_alphaii_Intercept
  alpha   <-matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2) 
  alpha   <-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
  
}
return(alpha)  
}



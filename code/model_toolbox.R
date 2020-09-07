#toolbox for example BEV competition
require(brms)
source("code/feasibility_toolbox.R")
source("code/growth_rates.R")





#function to generate a posterior growth rate (and an environmental growth rate) for each point of the posterior, for one model. exp_param is binary that tells it if to take into consideration if the model has an exponential parameter
posterior_parameters<-function(model, growth_fun, equilibrium_fun ,s ,g, exp_param){
  post        <- posterior_samples(model)
  
  lambda      <- exp(post$b_lambda_Intercept)
  lambda_env  <- exp(post$b_lambda_Intercept + post$b_lambda_env)
  
  alphaii <- post$b_alphaii_Intercept
  alphaii_env <-post$b_alphaii_Intercept + post$b_alphaii_env
  
  alphaij     <- post$b_alphaij_Intercept
  alphaij_env  <- post$b_alphaij_Intercept + post$b_alphaij_env
  
  if(exp_param){
    b           <- post$b_beta_Intercept
    b_env       <- post$b_beta_Intercept + post$b_beta_env
    growth      <- growth_fun(s,g,lambda,b)
    env_growth  <- growth_fun(s,g,lambda_env,b_env)
    equilibrium <- equilibrium_fun(s,g,lambda,b,alphaii)
    env_equilibrium <-equilibrium_fun(s,g,lambda_env,b_env,alphaii_env)
  }else{
    growth     <- growth_fun(s,g,lambda)
    env_growth <- growth_fun(s,g,lambda_env)
    equilibrium <- equilibrium_fun(s,g,lambda,alphaii)
    env_equilibrium <-equilibrium_fun(s,g,lambda_env,alphaii_env)
  }
  
 posterior <- as.data.frame(cbind(lambda, lambda_env, alphaii, alphaii_env, alphaij, alphaij_env,growth,env_growth, equilibrium, env_equilibrium))
 
 
  return(posterior)
}



#function to generate an alpha matrix based on one row of the posterior for each species. env is a binary that tells it if to take into consideration the environmental values of alpha
alpha_matrix <- function(vero_row, trcy_row, gi,gj, env){
  # vero_row <- as.list(vero_row)
  # trcy_row <- as.list(trcy_row)
  
  #if env, then alphas are alpha + alpha_env 
  if(env){
    alpha11 <- vero_row$b_alphaii_Intercept + vero_row$b_alphaii_env
    alpha21 <- trcy_row$b_alphaij_Intercept + trcy_row$b_alphaij_env
    alpha12 <- vero_row$b_alphaij_Intercept + vero_row$b_alphaij_env
    alpha22 <- trcy_row$b_alphaii_Intercept + trcy_row$b_alphaii_env
    alpha   <- matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2)
    alpha   <- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }else{
    alpha11 <- vero_row$b_alphaii_Intercept
    alpha21 <- trcy_row$b_alphaij_Intercept
    alpha12 <- vero_row$b_alphaij_Intercept
    alpha22 <- trcy_row$b_alphaii_Intercept
    alpha   <- matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2) 
    alpha   <- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }
  return(alpha)  
}


#function that spits out an omega, feasibility and theta calculated from the posterior of two models. env is a binary that tells it if to take iinto consideration the environmental variables
posterior_feasibility<-function(vero_post,trcy_post,si,sj,gi,gj,env){
  
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(num_posterior){
    
    omega_results       <-c()
    feasibility_results <-c()
    growth_results      <-c()
    theta_results       <-c()
    
    for( i in 1:nrow(vero_post)){
      
      #we get the corresponding posterior values, vero first, trcy second, gi (vero), gj(trcy)
      alpha  <- alpha_matrix(vero_post[i,],trcy_post[i,],gi,gj,env)
 
      if(env){
        r1 <- vero_post$env_growth[i]
        r2 <- trcy_post$env_growth[i]
      }else{
        r1 <- vero_post$growth[i]
        r2 <- trcy_post$growth[i]
      }
      
 
      #And estimate the feasability domain
      omega       <-Omega(alpha)
      feasibility <-test_feasibility(alpha,c(r1,r2))
      theta       <-theta(alpha,c(r1,r2))
      
      # #we save it 
      omega_results       <-c(omega_results,omega)
      feasibility_results <-c(feasibility_results, feasibility)
      theta_results       <-c(theta_results,theta)
    }
    
    
    pp<- cbind(omega_results,feasibility_results,theta_results)
    
    pp<-as.data.frame(pp)
    return(pp)
  }else{warning("Posterior distributions are not the same length")}
  
  
}

## function that warps the previous function to spit a data frame that has the posterior parameter values, as well as the omega values, theta values and feasibility for each row of the posterior. 
feasibility_wrapper<-function(vero_model,vero_function, vero_exp,trcy_model,trcy_function,trcy_exp,env){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  
  vero_post<- posterior_parameters(vero_model, fun = vero_function, si, gi, exp_param = vero_exp)
  trcy_post<- posterior_parameters(trcy_model, fun= trcy_function, sj,  gj, exp_param = trcy_exp)
  
  post <- posterior_feasibility(vero_post, trcy_post,si,sj,gi,gj, env = env)
  
  return(post)
}




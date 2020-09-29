#toolbox for example BEV competition
require(brms)
source("code/feasibility_toolbox.R")






#function to generate a posterior growth rate (and an environmental growth rate) for each point of the posterior, for one model. exp_param is binary that tells it if to take into consideration if the model has an exponential parameter
posterior_parameters <- function(model,s,g) {

    # extract samples from the model posterior for all parameters

    post        <- posterior_samples(model)

    # the lambda parameter is in log space and needs transformation
    lambda      <- exp(post$b_lambda_Intercept)
    lambda_env  <- exp(post$b_lambda_Intercept + post$b_lambda_env)

    # intraspecific alphas in control and treatment conditions
    alphaii     <- post$b_alphaii_Intercept
    alphaii_env <- post$b_alphaii_Intercept + post$b_alphaii_env

    # interspecific alphas in control and treatment conditions
    alphaij     <- post$b_alphaij_Intercept
    alphaij_env <- post$b_alphaij_Intercept + post$b_alphaij_env
    
    # some models and hence "growth rates" include additional parameters
    if("beta" %in% names(model$formula$pforms)) {
      beta     <- exp(post$b_beta_Intercept)
      beta_env <- exp(post$b_beta_Intercept + post$b_beta_env)
    }else{
      beta     <- NA
      beta_env <- NA
    }
    
    #These definitions are true for all models
    a <- lambda*g
    a_env <- lambda_env*g
    b <- 1 - ((1-g)*s)
    
    #to calculate different growth rates based on the model we define growth based on the name of the model
    if (model$name=="Beverton-Holt"){
      growth <- (a/b) - 1
      env_growth <- (a_env/b) - 1
    }
    
    if(model$name=="Lotka-Volterra"){
      growth <- 1 - (b/a) 
      env_growth <- 1 - (b/a_env)
    }
    
    if(model$name =="Ricker"){
      growth <- log(a/b)
      env_growth <- log(a_env/b)
    }
    
    if(model$name=="Hassell"){
      growth <- -1 + ((a/b)^(1/b))
      env_growth <-   -1 + ((a_env/b)^(1/b_env))
    }
    
    
    # calculate single-species equilibrium densities in control and treatment conditions
    equilibrium <- growth/(alphaii*g)
    env_equilibrium <- env_growth/(alphaii_env*g)
    
    # create a posterior sample of "more interpretable" model parameters
    posterior <-
      as.data.frame(
        cbind(
          lambda,
          lambda_env,
          alphaii,
          alphaii_env,
          alphaij,
          alphaij_env,
          growth,
          env_growth,
          equilibrium,
          env_equilibrium
        )
      )
    
         
    
     #because to rbind things they need to have the same number of columsn
      posterior$b     <- beta
      posterior$b_env <- beta_env
    
    
    # add in the other alpha parameter if it exists (for completeness)
    if("alphaik" %in% names(model$formula$pforms)) {
      posterior$alphaik     <- post$b_alphaii_Intercept
      posterior$alphaik_env <- post$b_alphaii_Intercept + post$b_alphaii_env
    }
    
    # add the model name to always be sure 
    posterior$model <- model$name
    return(posterior)
}


#function to generate an alpha matrix based on one row of the posterior for each species. env is a binary that tells it if to take into consideration the environmental values of alpha
alpha_matrix <- function(vero_row, trcy_row, gi,gj, env){
  # vero_row <- as.list(vero_row)
  # trcy_row <- as.list(trcy_row)
  
  #if env, then alphas are alpha + alpha_env 
  if(env){
    alpha11 <- vero_row$alphaii_env
    alpha21 <- trcy_row$alphaij_env
    alpha12 <- vero_row$alphaij_env
    alpha22 <- trcy_row$alphaii_env
    alpha   <- matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2) 
    alpha   <- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }else{
    alpha11 <- vero_row$alphaii
    alpha21 <- trcy_row$alphaij
    alpha12 <- vero_row$alphaij
    alpha22 <- trcy_row$alphaii
    alpha   <- matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2) 
    alpha   <- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }
  return(alpha)  
}


#function that spits out an omega, feasibility and theta calculated from the posterior of two models. env is a binary that tells it if to take iinto consideration the environmental variables
posterior_feasibility<-function(vero_post,trcy_post,gi,gj,env){
  
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(num_posterior){
    
    omega_results       <-c()
    feasibility_results <-c()
    growth_results      <-c()
    theta_results       <-c()
    
    for( i in 1:nrow(vero_post)){
      
      #we get the corresponding posterior values, vero first, trcy second, gi (vero), gj(trcy)
      alpha  <- alpha_matrix(vero_row=vero_post[i,],trcy_row =trcy_post[i,],gi=gi,gj=gj,env=env)
 
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
feasibility_wrapper<-function(vero_model,vero_function, vero_exp,vero_name,trcy_model,trcy_function,trcy_exp,trcy_name,env){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  
  vero_post<- posterior_parameters(model = vero_model, growth_fun = vero_function,s = si,g = gi,exp_param = vero_exp)
  trcy_post<- posterior_parameters(model = trcy_model, growth_fun = trcy_function,s = sj,g  = gj, exp_param = trcy_exp) 
  
  post <- posterior_feasibility(vero_post = vero_post,trcy_post = trcy_post,gi = gi,gj = gj,env = env)
  post$vero_model <- vero_name
  post$trcy_model <- trcy_name
  return(post)
}






# function to extract the main effects of a model fit
fixed_model<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  return(params)
}

#function to get the growth rate of a model, based on the mean parameter values. The model fit should have a name within, hopefully you have that when you read the model in a scipt!
get_fixed_growth<- function(model,s,g,env){
 
   params <- fixed_model(model)
  lambda <- exp(params$lambda_Intercept)
  lambda_env <- exp(params$lambda_Intercept + params$lambda_env)
  
  a <- lambda*g
  a_env <- lambda_env*g
  b <- 1 - ((1-g)*s)
  if("beta" %in% names(model$formula$pforms)) {
    beta     <- exp(params$beta_Intercept)
    beta_env <- exp(params$beta_Intercept + params$beta_env)
  }else{
    beta     <- NA
    beta_env <- NA
  }
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
     growth <- -1 + ((a/b)^(1/beta))
     env_growth <-   -1 + ((a_env/b)^(1/beta_env))
   }
  
  if(env){
    return(env_growth)
  } else{
    return(growth)
  } 
    }

#function to extract thealpha matrix based on two models, for the mean parameter values, env is a binary to tell it to take into consideration parameter values assosiated with the woody environment to calculate growth
get_fixed_alphas<-function(vero_model,gi,trcy_model,gj, env){
  vero<-fixed_model(vero_model)
  trcy<-fixed_model(trcy_model)
  
  if(env){
    
    a11<-vero$alphaii_Intercept + vero$alphaii_env
    a21<-trcy$alphaij_Intercept + trcy$alphaij_env
    a12<-vero$alphaij_Intercept + vero$alphaij_env
    a22<-trcy$alphaii_Intercept + trcy$alphaii_env
    
    alpha<-matrix( c(a11,a21,a12,a22) ,ncol=2,nrow=2) 
    alpha<- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
  }else{
    a11<-vero$alphaii_Intercept
    a21<-trcy$alphaij_Intercept
    a12<-vero$alphaij_Intercept
    a22<-trcy$alphaii_Intercept
    
    alpha<-matrix( c(a11,a21,a12,a22) ,ncol=2,nrow=2) 
    alpha<- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*") 
  }
  
  return(alpha)
}


#function to generate a posterior growth rate (and an environmental growth rate) for each point of the posterior, for one model. s and g are the seed survival and germination rates. Model is the brms model fit. It spits out the posterior parameters, the posterior growth and environmental growth and equilibrium abundances and environmental equilibrium abundances
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
      growth <- -1 + ((a/b)^(1/beta))
      env_growth <-   -1 + ((a_env/b)^(1/beta_env))
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
    
    
    
    #because to rbind things later on  they need to have the same number of columsn
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
posterior_feasibility<-function(vero_model,trcy_model,si,gi,sj,gj,env){
  
  #for the mean omega and theta we get the  alpha matrix for mean parameter values
  #for an environmental condition
  mean_alpha_matrix <- get_fixed_alphas(vero_model = vero_model,
                                        trcy_model = trcy_model,
                                        gi = gi,
                                        gj = gj,
                                        env = env)
  #as well as r1
  vero_growth <- get_fixed_growth(model = vero_model,
                                  s =si,
                                  g =gi, 
                                  env = env)
  
  #and r2
  trcy_growth <- get_fixed_growth( model = trcy_model,
                                   s = sj,
                                   g = gj,
                                   env =env)
  
  #and now we can calculate the mean
  
  Omega_mean <- Omega(mean_alpha_matrix)
  Theta_mean <- theta(mean_alpha_matrix, c(vero_growth, trcy_growth))
  
  # for the posterior feasibility
  trcy_post<-posterior_parameters(model = trcy_model, s = sj,g = gj)
  vero_post<-posterior_parameters(model = vero_model, s = si, g = gi)
  
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
    
    
    pp<- as.data.frame(cbind(omega_results,feasibility_results,theta_results))
    pp$Omega_mean <- Omega_mean
    pp$Theta_mean <- Theta_mean
    
    pp$vero_model <- vero_model$name
    pp$trcy_model <- trcy_model$name
    return(pp)
  }else{warning("Posterior distributions are not the same length")}
  
  
}





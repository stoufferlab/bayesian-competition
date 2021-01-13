
#Functions to extract parameters from bayesian fits and apply the feasibility framework 
#Most of these functions get used in the function posterior_feasibility

# function to extract the main effects of a model fit
fixed_model<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  return(params)
}

#function to get the growth rate of species given a model, based on the mean parameter values. The model (model) fit should have a name within, hopefully you have that when you read the model in a scipt! Env is a env is a binary to tell it to take into consideration parameter values assosiated with the woody environment to calculate growth.
get_fixed_growth<- function(model,
                            s,
                            g,
                            env){
 
  params <- fixed_model(model)
  lambda <- exp(params$lambda_Intercept)
  lambda_env <- exp(params$lambda_Intercept + params$lambda_env)
  
  a <- lambda*g
  a_env <- lambda_env*g
  b <- 1 - ((1-g)*s)
  
  #dealing with the hassell model
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

#function to extract the alpha matrix based on two models, for the mean parameter values, env is a binary to tell it to take into consideration parameter values assosiated with the woody environment to calculate growth
get_fixed_alphas<-function(vero_model,
                           gi,
                           trcy_model,
                           gj, 
                           env){
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
posterior_parameters <- function(model,
                                 s,
                                 g) {

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
alpha_matrix <- function(vero_row,
                         trcy_row,
                         gi,
                         gj, 
                         env){
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


#function that spits out an omega, feasibility and theta calculated from the posterior of two models. env is a binary that tells it if to take iinto consideration the environmental variables. Ni_max and Nj_max are the maximum expected abundances of species i (vero) and j (trcy)
posterior_feasibility <- function(vero_model,
                                  trcy_model,
                                  si,
                                  gi,
                                  Ni_max,
                                  sj,
                                  gj,
                                  Nj_max,
                                  env = FALSE,
                                  make_plot = FALSE){
  #for the mean omega and theta we get the  alpha matrix for mean parameter values
  #for an environmental condition either 0 (control) or 1 (woody)
  name <- paste0(vero_model$name, "and", trcy_model$name)
  print(name)
  mean_alpha_matrix <- get_fixed_alphas(
    vero_model = vero_model,
    trcy_model = trcy_model,
    gi = gi,
    gj = gj,
    env = env)
  

  
  #as well as r1 (vero's growth rate)
  vero_growth <- get_fixed_growth(
    model = vero_model,
    s = si,
    g = gi,
    env = env)
  
  #and r2 (trcy's growth rate)
  trcy_growth <- get_fixed_growth(
    model = trcy_model,
    s = sj,
    g = gj,
    env = env)
  #store them in a vector
  r <- c(vero_growth,
         trcy_growth)
  
  # Each model has its own constraints
  rconstraints <- list(
    lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
    upper = c(vero_model$constraints[2], trcy_model$constraints[2])
  )
  
  #And each species its maximum expected abundances
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  #Which determine the Radius for each alpha matrix used
  R <- determine_radius(alpha = mean_alpha_matrix,
                        Ni_max = Ni_max,
                        Nj_max = Nj_max)
  

  # Knowing these we can calculate the feasibility domain and its center for mean parameter values
  
    fixed_feasibility <- integrate_radii(
      alpha = mean_alpha_matrix,
      R = R,
      rconstraints = rconstraints,
      Nupper = Nupper)

  
  #Saaveda et al. estimation  
   fixed_feasibility_SA <- Omega_SA(alpha = mean_alpha_matrix)
   fixed_centroid_SA <- r_centroid(mean_alpha_matrix)
   fixed_theta_SA <-theta(r_c = fixed_centroid_SA,
                          r =r)
  # ou estimation of the center
   fixed_center <- r_feasible(
     alpha = mean_alpha_matrix,
     rconstraints = rconstraints,
     Nupper = Nupper,
     R_max = R,
     make_plot = make_plot)
   
  #check if our growth rates are feasible
    feasiblity_mean <- check_feasibility(
      r = r,
      alpha = mean_alpha_matrix,
      rconstraints = rconstraints,
      Nupper = Nupper )
  #Calculate the distance from the center
   distance_mean <- calculate_distance(center = fixed_center,
                                       r = r)
   theta_mean <- theta(r_c = fixed_center,
                       r = r )
  
  #we store the values of coexistence using the point estimates
  mean_parameters_results <- data.frame(
                                   "Omega_mean_saaveda"= fixed_feasibility_SA,
                                   "theta_mean_saavedra" = fixed_theta_SA,
                                   "Omega_mean"= fixed_feasibility,
                                    "distance_mean"=  distance_mean,
                                    "theta_mean" = theta_mean,
                                   "feasibility_mean"= feasiblity_mean,
                                   "R_mean"=R)
  print(mean_parameters_results)
  
  #######NOW for the posterior parameters#################################

  
  # we extract the posterior parameter values and growwth rates
  vero_post <- posterior_parameters(model = vero_model,
                                    s = si,
                                    g = gi)
  
  trcy_post <- posterior_parameters(model = trcy_model,
                                    s = sj,
                                    g = gj)
  
  
  # we check if they are the same size
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(!num_posterior){
    warning("Posterior distributions are not the same length")
  }else{
    
    print("working with the posterior distrubution")
    #just to work with them, should comment out this part aftewards
    vero_post<-vero_post[sample(nrow(vero_post), 1000), ]
    trcy_post<-trcy_post[sample(nrow(trcy_post), 1000), ]
    
    
    #to iterate over rows without using a loop
    x <- seq(1,nrow(vero_post),1) %>% as.list()
    
    posterior_parameters_results<-lapply(x,function(rows,gi,gj,rconstraints, Nupper){
      
      alpha  <- alpha_matrix(
          vero_row = vero_post[rows, ],
          trcy_row = trcy_post[rows, ],
          gi = gi,
          gj = gj,
          env = env )
      
      
      if (env) {
        r1 <- vero_post$env_growth[rows]
        r2 <- trcy_post$env_growth[rows]
        
      } else{
        r1 <- vero_post$growth[rows]
        r2 <- trcy_post$growth[rows]
      }
      
      r_post <- c(r1,r2)
      
     #we determine R for every alpha matrix
      R_post <- determine_radius(alpha = alpha, 
                                 Ni_max = Ni_max,
                                 Nj_max = Nj_max)
   
     
     
       omega_post <- integrate_radii(alpha = alpha,
                                     R = R_post,
                                     rconstraints = rconstraints,
                                     Nupper = Nupper )
    
      #Saavedras aproximation
      omega_post_SA <- Omega_SA(alpha = alpha)
      centroid_post_SA <- r_centroid(alpha=alpha)
      theta_post_SA <- theta(r_c = centroid_post_SA, 
                             r = r_post)
      #center of the domain
        center_post <- r_feasible(alpha = alpha,
                                  rconstraints = rconstraints,
                                  Nupper = Nupper,
                                  R_max = R_post,
                                  make_plot = FALSE)
      #are our growth rates feasible?
       feasibility_post <- check_feasibility(r= r_post,
                                             alpha = alpha,
                                             rconstraints = rconstraints,
                                             Nupper = Nupper )
      #how far away are they from the center
       distance_post <- calculate_distance(center = center_post,
                                           r = r_post)
       theta_post <- theta(r_c = center_post,
                           r = r_post)
      #all togethe
      post_results <- data.frame(
                                 "Omega_saaveda"= omega_post_SA,
                                 "theta_saavedra"= theta_post_SA,
                                 "Omega"= omega_post, 
                                  "distance"= distance_post,
                                 "theta" = theta_post,
                                 "feasibility" = feasibility_post,
                                 "Radius" = R_post)
   
       print(x)
      return(post_results)
      
    }, gi = gi,
    gj = gj,
    rconstraints = rconstraints, 
    Nupper = Nupper)
   
    posterior_parameters_results <- do.call(rbind, posterior_parameters_results)
    
    all_results <- cbind(mean_parameters_results, posterior_parameters_results)
    all_results[,"vero_model"] <- vero_model$name
    all_results[,"trcy_model"] <- trcy_model$name
   # print(all_results)
     file <- paste0(name,".RDS")
    # print(file)
    saveRDS(object = all_results,file = paste0(name,".RDS"))
    return(all_results) 
  }
    }

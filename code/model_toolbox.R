
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
                                  bounded = TRUE){
  #for the mean omega and theta we get the  alpha matrix for mean parameter values
  #for an environmental condition either 0 (control) or 1 (woody)
  name <- paste0(vero_model$name, "&", trcy_model$name)
  print(name)
  #we determine the alpha matrix for the mean parameter values
  alpha_mean <- get_fixed_alphas(
    vero_model = vero_model,
    trcy_model = trcy_model,
    gi = gi,
    gj = gj,
    env = env)
  #as well as r1 (vero's growth rate)
  vero_growth_mean <- get_fixed_growth(
    model = vero_model,
    s = si,
    g = gi,
    env = env)
  #and r2 (trcy's growth rate)
  trcy_growth_mean <- get_fixed_growth(
    model = trcy_model,
    s = sj,
    g = gj,
    env = env)
  #store them in a vector
  r_mean <- c(vero_growth_mean,
              trcy_growth_mean)
  print(r_mean)
  #determine if there are boundaries on the growth rates
  if(!bounded) {
    # Each model has its own constraints
    rconstraints <- list(
      lower = c(-Inf, -Inf),
      upper = c(Inf, Inf)
    )
    
  }else{
    # Each model has its own constraints
    rconstraints <- list(
      lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
      upper = c(vero_model$constraints[2], trcy_model$constraints[2])
    )  
  }

  
  
  #And each species its maximum expected abundances
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  #Which determine the Radius for each alpha matrix used
  R_mean <- determine_radius(alpha = alpha_mean,
                             Ni_max = Ni_max,
                             Nj_max = Nj_max)
  print(R_mean)
  
  
  #we plot how it looks just for kicks
  # png(file= paste(name,"_",env,"_",bounded,"_","mean", ".png"))
  # plot(0,0,
  #      xlim= c(-R_mean, R_mean),
  #      ylim= c(-R_mean, R_mean),
  #      # xlim=c(min(bounds_mean$ri), max(bounds_mean$ri)),
  #      # ylim=c(min(bounds_mean$rj), max(bounds_mean$rj)),
  #      type='n',
  #      xlab=expression(italic(r[i])),
  #      ylab=expression(italic(r[j]))
  # )
  # abline(h=0,lty='dashed',lwd=1.5)
  # abline(v=0,lty='dashed',lwd=1.5)
  #Saaveda et al. estimation, to compare
  omega_SA_mean <- Omega_SA(alpha = alpha_mean)
  
  centroid_SA_mean <- r_centroid(alpha_mean)
  
  theta_SA_mean <- theta(r_c = centroid_SA_mean,
                         r = r_mean)
  
  feasibility_SA_mean <- test_feasibility_saavedra(alpha = alpha_mean,
                                                   r = r_mean)
  #check if our growth rates are feasible, accordig to our constraints
  feasiblity_mean <- check_point(r =r_mean,
                                 R_max = R_mean,
                                 inv_alpha =  inverse_matrix(alpha_mean),
                                 rconstraints = rconstraints,
                                 Nupper = Nupper)
  feas_na <- is.na(feasiblity_mean)
  #in this case if it is NA it means it is unfeasible
  feasiblity_mean <-  ifelse(feas_na,0, feasiblity_mean)
  
  
  #Now we determine the feasibility shape with our Monte Carlo Sampling
 integration_mean<- integrate_area(R_max = R_mean,
                 alpha = alpha_mean,
                 rconstraints = rconstraints,
                 Nupper = Nupper,
                 desired_feasible = 1000,
                 max_samples = 1e6
  )
 #which spits out the propotion of the area that is feasible, or the feasibility domain
 Omega_mean <- integration_mean$proportion
 #And also the coordinates of all the points that are feasible
 shape_mean <- integration_mean$coords
 #and the points that are unfeasible
 unfeasible_mean <- integration_mean$unfeasible
 # which tell us the bounds of the feasibility domain
 shape_bounds_mean <- determine_boundary_shape(shape = shape_mean)
 bounds_mean <- shape_bounds_mean$bounds
 # and also the area the area of the feasibility domain
 area_mean <- shape_bounds_mean$area
 #with the bounds we can then get the distance from the limit of our growth rates
 distances_mean <- distance_from_limit(r=r_mean,
                                       shape = bounds_mean,
                                       feasibility = feasiblity_mean)
 distance_growth_mean <- distances_mean$growth_distance
 distance_center_mean <- distances_mean$center_distance
 #but also we get the proportion of things inside the convex hull
 convex_mean <- calculate_convex(shape = bounds_mean,
                                 unfeasible = unfeasible_mean)
 prop <-nrow(shape_mean)/nrow(unfeasible_mean)
 # And we keep track of everything
 #we store the values of coexistence using the point estimates
 mean_parameters_results <- data.frame(
   "Omega_saavedraa_mean"= omega_SA_mean,
   "theta_saavedra_mean" = theta_SA_mean,
   "feasibility_saavedra_mean" = feasibility_SA_mean,
   "Omega_mean"= Omega_mean,
   "area_mean"= area_mean,
   "distance_center_mean"=  distance_center_mean,
   "distance_growth_mean"=  distance_growth_mean,
   "feasibility_mean"= feasiblity_mean,
   "R_mean"=R_mean,
   "convex_mean"= convex_mean,
   "our_proportion"=prop)
 print(mean_parameters_results)
 
# 
#  col1 <- rethinking::col.alpha("mediumseagreen", alpha=0.1)
#  col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
#  lines(bounds_mean$ri, bounds_mean$rj, col= "black", lwd=2)
#  points(shape_mean$ri, shape_mean$rj, pch=20, col=col1)
#  points(unfeasible_mean$ri, unfeasible_mean$rj, pch=20, col=col2)
 # dev.off()

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
    vero_post_sample<-vero_post[1:200, ]
    trcy_post_sample<-trcy_post[1:200, ]
    
    #to iterate over rows without using a loop
    x <- seq(1,nrow(vero_post_sample),1) %>% as.list()

    posterior_parameters_results<-lapply(x,function(rows,
                                                    vero_post_sample,
                                                    trcy_post_sample,
                                                    gi,
                                                    gj,
                                                    rconstraints,
                                                    Nupper){
 print("WORKING WITH THE POSTEIOR")
      print(rows)
      #we get the alpha matrix for a point in the posteior and the growth rates
      alpha  <- alpha_matrix(
        vero_row = vero_post_sample[rows, ],
        trcy_row = trcy_post_sample[rows, ],
        gi = gi,
        gj = gj,
        env = env )
   print(alpha)   
      if (env) {
        r1 <- vero_post$env_growth[rows]
        r2 <- trcy_post$env_growth[rows]
        
      } else{
        r1 <- vero_post$growth[rows]
        r2 <- trcy_post$growth[rows]
      }
      
      r_post <- c(r1,r2)
      print("growth rates")
      print(r_post)
      #we determine R for every alpha matrix
      R_post <- determine_radius(alpha = alpha, 
                                 Ni_max = Ni_max,
                                 Nj_max = Nj_max)
      
      # png(file= paste(name,"_",env,"_",bounded,"_",rows, ".png"))
      # plot(0,0,
      #      xlim= c(-R_post, R_post),
      #      ylim= c(-R_post, R_post),
      #      # xlim=c(min(bounds_mean$ri), max(bounds_mean$ri)),
      #      # ylim=c(min(bounds_mean$rj), max(bounds_mean$rj)),
      #      type='n',
      #      xlab=expression(italic(r[i])),
      #      ylab=expression(italic(r[j]))
      # )
      # abline(h=0,lty='dashed',lwd=1.5)
      # abline(v=0,lty='dashed',lwd=1.5)
      
      
      
      print("Radius")
      print(R_post)
      #Saavedras aproximation
      omega_post_SA <- Omega_SA(alpha = alpha)
      
      centroid_post_SA <- r_centroid(alpha = alpha)
      
      theta_post_SA <- theta(r_c = centroid_post_SA,
                             r = r_post)
      
      feasibility_post_SA <- test_feasibility_saavedra(alpha = alpha, 
                                                       r = r_post)
      
      #our approximation ############################################
      f_post <- check_point(r =r_post,
                                     R_max = R_post,
                                     inv_alpha =  inverse_matrix(alpha),
                                     rconstraints = rconstraints,
                                     Nupper = Nupper)
      
      feas_post <- is.na(f_post)
      print(feas_post)
      #in this case if it is NA it means it is unfeasible
      f_post <-  ifelse(feas_post,0, f_post)
     
      #Now we determine the feasibility shape with our Monte Carlo Integration
      integration_post <- integrate_area(R_max = R_post,
                                        alpha = alpha,
                                        rconstraints = rconstraints,
                                        Nupper = Nupper,
                                        desired_feasible = 1000,
                                        max_samples = 1e6
      )
      
      #which spits out the propotion of the area that is feasible, or the feasibility domain
      Omega_post <- integration_post$proportion
      #And also the coordinates of all the points that are feasible
      shape_post <- integration_post$coords
      # which tell us the bounds of the feasibility domain
      shape_bounds_post <- determine_boundary_shape(shape = shape_post)
      bounds_post <- shape_bounds_post$bounds
      # and also the area the area of the feasibility domain
      area_post <- shape_bounds_post$area
      
      unfeasible_post <- integration_post$unfeasible
      #with the bounds we can then get the distance from the limit of our growth rates
      distances_post <- distance_from_limit(r=r_post,
                                            shape = bounds_post,
                                            feasibility = f_post)
      print("distances")
      
      convex_post <- calculate_convex(shape = bounds_post,
                                      unfeasible = unfeasible_post)
      
      print(distances_post)
      distance_growth_post <- distances_post$growth_distance
      distance_center_post <- distances_post$center_distance
     
      
      # 
      # col1 <- rethinking::col.alpha("mediumseagreen", alpha=0.1)
      # col2 <- rethinking::col.alpha("grey50", alpha = 0.1)
      # 
      # lines(bounds_post$ri, bounds_post$rj, col= "black", lwd=2)
      # points(shape_post$ri, shape_post$rj, pch=20, col=col1)
      # points(unfeasible_post$ri, unfeasible_post$rj, pch=20, col=col2)
      # dev.off()
      #all together, they make a row of results
      post_results <- data.frame(
        "Omega_saaveda"= omega_post_SA,
        "theta_saavedra"= theta_post_SA,
        "feasibility_saaveda"= feasibility_post_SA,
        "Omega"= Omega_post, 
        "area" = area_post,
        "distance_center"= distance_center_post,
        "distance_growth"= distance_growth_post,
        "feasibility" = f_post,
        "Radius" = R_post,
        "convex"= convex_post,
        "r1" = r1,
        "r2"= r2,
        "alpha11"= alpha[1,1],
        "alpha21"= alpha[2,1],
        "alpha12"= alpha[1,2],
        "alpha22"= alpha[2,2])
      

      return(post_results)
      
    }, vero_post_sample = vero_post_sample ,
    trcy_post_sample= trcy_post_sample,
    gi = gi,
    gj = gj,
    rconstraints = rconstraints, 
    Nupper = Nupper)
    
    posterior_parameters_results <- do.call(rbind, posterior_parameters_results)
    
    all_results <- cbind(mean_parameters_results, posterior_parameters_results)
    all_results[,"vero_model"] <- vero_model$name
    all_results[,"trcy_model"] <- trcy_model$name
    
   # file <- paste0(name,".RDS")
    # print(file)
   # saveRDS(object = all_results,file = paste0(name,".RDS")
    return(all_results) 
  }
}

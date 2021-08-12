
##function to show variation in rho and delta including the uncertainty of an individual parameter


posterior_feasibility_sensitivity <- function(vero_model,
                                  trcy_model,
                                  si,
                                  gi,
                                  Ni_max,
                                  sj,
                                  gj,
                                  Nj_max,
                                  env = FALSE){
  
  
  rconstraints <- list(
    lower = c(vero_model$constraints[1], trcy_model$constraints[1]),
    upper = c(vero_model$constraints[2], trcy_model$constraints[2])
  )  
  
  
  #And each species its maximum expected abundances
  Nupper <- c(i = Ni_max/gi,
              j = Nj_max/gj)
  
  name <- paste0(vero_model$name, "&", trcy_model$name)
  
  alpha_mean <- get_fixed_alphas(
    vero_model = vero_model,
    trcy_model = trcy_model,
    gi = gi,
    gj = gj,
    env = env)
  
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
  
  
  # we extract the posterior parameter values and growwth rates
  vero_post <- posterior_parameters(model = vero_model,
                                    s = si,
                                    g = gi) %>% select(lambda,
                                                       lambda_env,
                                                       alphaii,
                                                       alphaii_env,
                                                       alphaij,
                                                       alphaij_env)
  
  trcy_post <- posterior_parameters(model = trcy_model,
                                    s = sj,
                                    g = gj)%>% select(lambda,
                                                      lambda_env,
                                                      alphaii,
                                                      alphaii_env,
                                                      alphaij,
                                                      alphaij_env)
  
  vero_post <- vero_post[1:500,]
  trcy_post <- trcy_post[1:500,]
  
  
  if(!env){
    parameters <- c("lambda",
                    "alphaii",
                    "alphaij")
  #vero  
  parameter_sensitivity_vero <- lapply(parameters, function(param){
    
      vero_name <-vero_model$name
      varying_parameter <- select(vero_post, eval(param) )
      
      results <-c()
      
      if(grepl("lambda", param)){

        for(i in 1:nrow(varying_parameter)){
          
          lambda_i <- varying_parameter[i,]
          a <- lambda_i*gi
          b <- 1 - ((1-gi)*si)
          
          if (vero_name=="Beverton-Holt"){
            growth <- (a/b) - 1
          }
          
          if(vero_name =="Ricker"){
            growth <- log(a/b)
          }
          
          R <- determine_radius(alpha = alpha_mean,
                                Ni_max = Ni_max, 
                                Nj_max = Nj_max)
          r <- c(growth, trcy_growth_mean)
          
          results_row<-  structural_stability_wrapper(R= R,
                                       alpha = alpha_mean,
                                       rconstraints = rconstraints,
                                       Nupper = Nupper,
                                       r =r )
          
          
          results <- rbind(results, results_row)
            
        }#end of i 
       
      }#end of lambda
      
      if(grepl("alpha", param)){
        #check out inf it is intra or inter
       which_alpha<- str_split_fixed(param, "i", 2)[2]
  
       #if this is intraspecific
       if(which_alpha=="i"){
         #then this is alpha11
         for(i in 1:nrow(varying_parameter)){
           
           alpha_mean[1,1] <- varying_parameter[i,]
           
           
           R <- determine_radius(alpha = alpha_mean,
                                 Ni_max = Ni_max, 
                                 Nj_max = Nj_max)
           r <- c(vero_growth_mean, trcy_growth_mean)
           
           results_row<-  structural_stability_wrapper(R= R,
                                                       alpha = alpha_mean,
                                                       rconstraints = rconstraints,
                                                       Nupper = Nupper,
                                                       r =r )
           
           
           results <- rbind(results, results_row)
           
         }#end of i
       } # end of intraspecific
        
       if(which_alpha=="j"){
         #then this is alpha12
         for(i in 1:nrow(varying_parameter)){
           
           alpha_mean[1,2] <- varying_parameter[i,]
           
           
           R <- determine_radius(alpha = alpha_mean,
                                 Ni_max = Ni_max, 
                                 Nj_max = Nj_max)
           r <- c(vero_growth_mean, trcy_growth_mean)
           
           results_row<-  structural_stability_wrapper(R= R,
                                                       alpha = alpha_mean,
                                                       rconstraints = rconstraints,
                                                       Nupper = Nupper,
                                                       r =r )
           
           
           results <- rbind(results, results_row)
           
         }#end of i
       } # end of interespecific
       
        
      }#end of alpha
      
      
      if(param == "lambda"){
        results$parameter <- "lambda_i"
      }
      if(param == "alphaii"){
        results$parameter <- "alphaii"
      }
      if(param == "alphaij"){
        results$parameter <- "alphaij"
      }
      
      return(results)
      
    })#end of vero
  parameter_sensitivity_vero <- do.call(rbind, parameter_sensitivity_vero)  
    
  #trcy  
  parameter_sensitivity_trcy <- lapply(parameters, function(param){
    
    trcy_name <-trcy_model$name
    varying_parameter <- select(trcy_post, eval(param) )
    
    results <-c()
    
    if(grepl("lambda", param)){
      
      for(i in 1:nrow(varying_parameter)){
        
        lambda_i <- varying_parameter[i,]
        a <- lambda_i*gi
        b <- 1 - ((1-gi)*si)
        
        if (trcy_name=="Beverton-Holt"){
          growth <- (a/b) - 1
        }
        
        if(trcy_name =="Ricker"){
          growth <- log(a/b)
        }
        
        R <- determine_radius(alpha = alpha_mean,
                              Ni_max = Ni_max, 
                              Nj_max = Nj_max)
        r <- c(vero_growth_mean, growth)
        
        results_row<-  structural_stability_wrapper(R= R,
                                                    alpha = alpha_mean,
                                                    rconstraints = rconstraints,
                                                    Nupper = Nupper,
                                                    r =r )
        
        
        results <- rbind(results, results_row)
        
      }#end of i 
      
    }#end of lambda
    
    if(grepl("alpha", param)){
      #check out inf it is intra or inter
      which_alpha<- str_split_fixed(param, "i", 2)[2]
      
      #if this is intraspecific
      if(which_alpha=="i"){
        #then this is alpha22
        for(i in 1:nrow(varying_parameter)){
          
          alpha_mean[2,2] <- varying_parameter[i,]
          
          
          R <- determine_radius(alpha = alpha_mean,
                                Ni_max = Ni_max, 
                                Nj_max = Nj_max)
          r <- c(vero_growth_mean, trcy_growth_mean)
          
          results_row<-  structural_stability_wrapper(R= R,
                                                      alpha = alpha_mean,
                                                      rconstraints = rconstraints,
                                                      Nupper = Nupper,
                                                      r =r )
          
          
          results <- rbind(results, results_row)
          
        }#end of i
      } # end of intraspecific
      
      if(which_alpha=="j"){
        #then this is alpha12
        for(i in 1:nrow(varying_parameter)){
          
          alpha_mean[2,1] <- varying_parameter[i,]
          
          
          R <- determine_radius(alpha = alpha_mean,
                                Ni_max = Ni_max, 
                                Nj_max = Nj_max)
          r <- c(vero_growth_mean, trcy_growth_mean)
          
          results_row<-  structural_stability_wrapper(R= R,
                                                      alpha = alpha_mean,
                                                      rconstraints = rconstraints,
                                                      Nupper = Nupper,
                                                      r =r )
          
          
          results <- rbind(results, results_row)
          
        }#end of i
      } # end of interespecific
      
      
    }#end of alpha
    
    if(param == "lambda"){
      results$parameter <- "lambda_j"
    }
    if(param == "alphaii"){
      results$parameter <- "alphajj"
    }
    if(param == "alphaij"){
      results$parameter <- "alphaji"
    }
    
    return(results)
    
  })#end of trcy
  parameter_sensitivity_trcy <- do.call(rbind, parameter_sensitivity_trcy)  
  
  parameter_sensitivity <- rbind(parameter_sensitivity_vero, parameter_sensitivity_trcy)
  
  parameter_sensitivity$vero_model <- vero_model$name
  parameter_sensitivity$trcy_model <- trcy_model$name
  parameter_sensitivity$environment <- "open"
    
  }else{
    parameters <- c("lambda_env",
                    "alphaii_env",
                    "alphaij_env")
    #vero  
    parameter_sensitivity_vero <- lapply(parameters, function(param){
      
      vero_name <-vero_model$name
      varying_parameter <- select(vero_post, eval(param) )
      
      results <-c()
      
      if(grepl("lambda", param)){
        
        for(i in 1:nrow(varying_parameter)){
          
          lambda_i <- varying_parameter[i,]
          a <- lambda_i*gi
          b <- 1 - ((1-gi)*si)
          
          if (vero_name=="Beverton-Holt"){
            growth <- (a/b) - 1
          }
          
          if(vero_name =="Ricker"){
            growth <- log(a/b)
          }
          
          R <- determine_radius(alpha = alpha_mean,
                                Ni_max = Ni_max, 
                                Nj_max = Nj_max)
          r <- c(growth, trcy_growth_mean)
          
          results_row<-  structural_stability_wrapper(R= R,
                                                      alpha = alpha_mean,
                                                      rconstraints = rconstraints,
                                                      Nupper = Nupper,
                                                      r =r )
          
          
          results <- rbind(results, results_row)
          
        }#end of i 
        
      }#end of lambda
      
      if(grepl("alpha", param)){

        #check out inf it is intra or inter
        which_alpha<- str_split_fixed(param, "i", 2)[2]
        
        #if this is intraspecific
        if(which_alpha=="i_env"){
          #then this is alpha11
          for(i in 1:nrow(varying_parameter)){
            
            alpha_mean[1,1] <- varying_parameter[i,]
            
            
            R <- determine_radius(alpha = alpha_mean,
                                  Ni_max = Ni_max, 
                                  Nj_max = Nj_max)
            r <- c(vero_growth_mean, trcy_growth_mean)
            
            results_row<-  structural_stability_wrapper(R= R,
                                                        alpha = alpha_mean,
                                                        rconstraints = rconstraints,
                                                        Nupper = Nupper,
                                                        r =r )
            
            
            results <- rbind(results, results_row)
            
          }#end of i
        } # end of intraspecific
        
        if(which_alpha=="j_env"){
          #then this is alpha12
          for(i in 1:nrow(varying_parameter)){
            
            alpha_mean[1,2] <- varying_parameter[i,]
            
            
            R <- determine_radius(alpha = alpha_mean,
                                  Ni_max = Ni_max, 
                                  Nj_max = Nj_max)
            r <- c(vero_growth_mean, trcy_growth_mean)
            
            results_row<-  structural_stability_wrapper(R= R,
                                                        alpha = alpha_mean,
                                                        rconstraints = rconstraints,
                                                        Nupper = Nupper,
                                                        r =r )
            
            
            results <- rbind(results, results_row)
            
          }#end of i
        } # end of interespecific
        
        
      }#end of alpha
      
      
      if(param == "lambda_env"){
        results$parameter <- "lambda_i"
      }
      if(param == "alphaii_env"){
        results$parameter <- "alphaii"
      }
      if(param == "alphaij_env"){
        results$parameter <- "alphaij"
      }
      
      return(results)
      
    })#end of vero
    parameter_sensitivity_vero <- do.call(rbind, parameter_sensitivity_vero)  
    
    #trcy  
    parameter_sensitivity_trcy <- lapply(parameters, function(param){
      
      trcy_name <-trcy_model$name
      varying_parameter <- select(trcy_post, eval(param) )
      
      results <-c()
      
      if(grepl("lambda", param)){
        
        for(i in 1:nrow(varying_parameter)){
          
          lambda_i <- varying_parameter[i,]
          a <- lambda_i*gi
          b <- 1 - ((1-gi)*si)
          
          if (trcy_name=="Beverton-Holt"){
            growth <- (a/b) - 1
          }
          
          if(trcy_name =="Ricker"){
            growth <- log(a/b)
          }
          
          R <- determine_radius(alpha = alpha_mean,
                                Ni_max = Ni_max, 
                                Nj_max = Nj_max)
          r <- c(vero_growth_mean, growth)
          
          results_row<-  structural_stability_wrapper(R= R,
                                                      alpha = alpha_mean,
                                                      rconstraints = rconstraints,
                                                      Nupper = Nupper,
                                                      r =r )
          
          
          results <- rbind(results, results_row)
          
        }#end of i 
        
      }#end of lambda
      
      if(grepl("alpha", param)){
        #check out inf it is intra or inter
        which_alpha<- str_split_fixed(param, "i", 2)[2]
        
        #if this is intraspecific
        if(which_alpha=="i_env"){
          #then this is alpha22
          for(i in 1:nrow(varying_parameter)){
            
            alpha_mean[2,2] <- varying_parameter[i,]
            
            
            R <- determine_radius(alpha = alpha_mean,
                                  Ni_max = Ni_max, 
                                  Nj_max = Nj_max)
            r <- c(vero_growth_mean, trcy_growth_mean)
            
            results_row<-  structural_stability_wrapper(R= R,
                                                        alpha = alpha_mean,
                                                        rconstraints = rconstraints,
                                                        Nupper = Nupper,
                                                        r =r )
            
            
            results <- rbind(results, results_row)
            
          }#end of i
        } # end of intraspecific
        
        if(which_alpha=="j_env"){
          #then this is alpha12
          for(i in 1:nrow(varying_parameter)){
            
            alpha_mean[2,1] <- varying_parameter[i,]
            
            
            R <- determine_radius(alpha = alpha_mean,
                                  Ni_max = Ni_max, 
                                  Nj_max = Nj_max)
            r <- c(vero_growth_mean, trcy_growth_mean)
            
            results_row<-  structural_stability_wrapper(R= R,
                                                        alpha = alpha_mean,
                                                        rconstraints = rconstraints,
                                                        Nupper = Nupper,
                                                        r =r )
            
            
            results <- rbind(results, results_row)
            
          }#end of i
        } # end of interespecific
        
        
      }#end of alpha
      
      if(param == "lambda_env"){
        results$parameter <- "lambda_j"
      }
      if(param == "alphaii_env"){
        results$parameter <- "alphajj"
      }
      if(param == "alphaij_env"){
        results$parameter <- "alphaji"
      }
      
      return(results)
      
    })#end of trcy
    parameter_sensitivity_trcy <- do.call(rbind, parameter_sensitivity_trcy)  
    
    parameter_sensitivity <- rbind(parameter_sensitivity_vero, parameter_sensitivity_trcy)
    
    parameter_sensitivity$vero_model <- vero_model$name
    parameter_sensitivity$trcy_model <- trcy_model$name
    parameter_sensitivity$environment <- "woody"
    
  }
 
  return(parameter_sensitivity)

}
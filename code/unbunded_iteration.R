#iterating over models

library(brms)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#We source everything known to human kind...


source("code/read_models.R")
source("code/model_toolbox.R")
source("code/model_combo.R")
source("code/integration_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")
source("code/determine_boundary.R")


#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033



posterior_feasibility_unbounded <- function(vero_model,
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
  
  # Each model has its own constraints. IN THIS CASE THEY ARE UNBOUNDED
  rconstraints <- list(
    lower = c(-Inf, -Inf),
    upper = c(Inf, Inf)
  )
  
  #And each species its maximum expected abundances
  Nupper <- c(i = Ni_max,
              j = Nj_max)
  #int
  R <- 1
  
  
  # Knowing these we can calculate the feasibility domain and its center for mean parameter values
  
  fixed_omega <- integrate_radii(
    alpha = mean_alpha_matrix,
    R = R,
    rconstraints = rconstraints,
    Nupper = Nupper)
  
  
  #Saaveda et al. estimation, to compare
  fixed_omega_SA <- Omega_SA(alpha = mean_alpha_matrix)
  fixed_centroid_SA <- r_centroid(mean_alpha_matrix)
  fixed_theta_SA <- theta(r_c = fixed_centroid_SA,
                          r = r)
  fixed_feasibility_SA <- test_feasibility_saavedra(alpha = mean_alpha_matrix,
                                                    r = r)
  #check if our growth rates are feasible
  fixed_feasiblity <- check_feasibility(
    r = r,
    alpha = mean_alpha_matrix,
    rconstraints = rconstraints,
    Nupper = Nupper )
  
  
  # ou estimation of the distance from the bounds
  fixed_distance<-  distance_from_limit(alpha = mean_alpha_matrix,
                                        R_max = R,
                                        rconstraints = rconstraints, 
                                        Nupper = Nupper,
                                        r = r,
                                        feasibility = fixed_feasiblity)
  
  
  
  
  #we store the values of coexistence using the point estimates
  mean_parameters_results <- data.frame(
    "Omega_mean_saaveda"= fixed_omega_SA,
    "theta_mean_saavedra" = fixed_theta_SA,
    "feasibility_mean_saavedra" = fixed_feasibility_SA,
    "Omega_mean"= fixed_omega,
    "distance_mean"=  fixed_distance,
    "feasibility_mean"= fixed_feasiblity,
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
    vero_post<-vero_post[sample(nrow(vero_post), 10), ]
    trcy_post<-trcy_post[sample(nrow(trcy_post), 10), ]
    
    
    #to iterate over rows without using a loop
    x <- seq(1,nrow(vero_post),1) %>% as.list()
    
    posterior_parameters_results<-lapply(x,function(rows,gi,gj,rconstraints, Nupper){
      
      alpha  <- alpha_matrix(
        vero_row = vero_post[rows, ],
        trcy_row = trcy_post[rows, ],
        gi = gi,
        gj = gj,
        env = env )
      print("this particula alpha matrix is")
      print(alpha)
      if (env) {
        r1 <- vero_post$env_growth[rows]
        r2 <- trcy_post$env_growth[rows]
        
      } else{
        r1 <- vero_post$growth[rows]
        r2 <- trcy_post$growth[rows]
      }
      
      r_post <- c(r1,r2)
      print("the growth rate vectir is")
      print(r_post)
      #we determine R for every alpha matrix
      R_post <- 1
      print("Radius is")
      print(R_post)
      #and the size of the feasibility domain
      omega_post <- integrate_radii(alpha = alpha,
                                    R = R_post,
                                    rconstraints = rconstraints,
                                    Nupper = Nupper )
      print("the feasibility domain is") 
      print(omega_post)
      
      #are our growth rates feasible?
      feasibility_post <- check_feasibility(r= r_post,
                                            alpha = alpha,
                                            rconstraints = rconstraints,
                                            Nupper = Nupper )
      
      print("distance is")
      distance_post <-  distance_from_limit(alpha = alpha,
                                            R_max = R_post,
                                            rconstraints = rconstraints, 
                                            Nupper = Nupper,
                                            r = r_post,
                                            feasibility = feasibility_post)
      print(distance_post)
      #Saavedras aproximation
      omega_post_SA <- Omega_SA(alpha = alpha)
      centroid_post_SA <- r_centroid(alpha = alpha)
      theta_post_SA <- theta(r_c = centroid_post_SA,
                             r = r_post)
      feasibility_post_SA <- test_feasibility_saavedra(alpha = alpha, 
                                                       r = r_post)
      
      #all togethe
      post_results <- data.frame(
        "Omega_saaveda"= omega_post_SA,
        "theta_saavedra"= theta_post_SA,
        "feasibility_saaveda"= feasibility_post_SA,
        "Omega"= omega_post, 
        "distance"= distance_post,
        "feasibility" = feasibility_post,
        "Radius" = R_post)
      
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
    # saveRDS(object = all_results,file = paste0(name,".RDS"))
    return(all_results) 
  }
}




fixed_row_model <- function(vero_models,
                            trcy_model,
                            si=si,
                            gi = gi,
                            Ni = Ni,
                            sj = sj,
                            gj = gj,
                            Nj = Nj,
                            make_plot=FALSE,
                            env = FALSE){
  one_row<-lapply(vero_models,function(m){
    post <- posterior_feasibility_unbounded(vero_model=m,
                                  trcy_model = trcy_model,
                                  si = si,
                                  gi = gi,
                                  sj = sj,
                                  gj = gj,
                                  Ni = Ni,
                                  Nj = Nj,
                                  env = env,
                                  make_plot = make_plot)
    
    return(post)
  })
  
  
  all_posteriors <- do.call(rbind, one_row)
  return(all_posteriors)
  
}

#trcy models is a list of models as well


combined_models<-function(vero_models,
                          trcy_models,
                          si=si,
                          gi = gi,
                          Ni = Ni,
                          sj = sj,
                          gj = gj,
                          Nj = Nj,
                          make_plot=FALSE,
                          env = FALSE){
  many_rows <- lapply(trcy_models, function(m2){
    one_row<-fixed_row_model(vero_models = vero_models,
                             trcy_model = m2,
                             si=si,
                             gi = gi,
                             Ni = Ni,
                             sj = sj,
                             gj = gj,
                             Nj = Nj,
                             make_plot=make_plot,
                             env = env)
    return(one_row)
  })
  
  all_rows<- do.call(rbind, many_rows)
  return(all_rows)
}



vero_models <- list( vero_bh_multispecies_poisson.rds,
                     vero_lv_multispecies_poisson.rds,
                     vero_rc_multispecies_poisson.rds)

trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                   trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds)

model_grid_sunny<- combined_models(vero_models = vero_models,
                                   trcy_models = trcy_models,
                                   si =si,
                                   gi =gi,
                                   gj =gj,
                                   sj=sj,
                                   Ni = 1e4,
                                   Nj =1e4,
                                   env=FALSE,
                                   make_plot = FALSE)

saveRDS(model_grid_sunny,
        file = "results_sunny_unbounded_22jan21.RDS")




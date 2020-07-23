#toolbox for example BEV competition
require(brms)
require(tidyverse)
require(ggplot2)
require(rethinking)
source("code/feasibility_toolbox.R")
source("code/extract_params.R")
source("code/growth_rates.R")
source("code/model_growth_rates.R")



#This function return a data frame with omegas, thetas, growth rates  based on the posterior distribution of parameter values given:
#vero_model is a model fit for vero (brms object)
#trcy_model is a model fit for trcy (brms object)
#vero_growth is a tibble that has the growth rates of vero for different models based on the posterior of  lambdai
#trcy_growth is a tibble that has the growth rates of trcy for different models based on the posterior of lambdaj
#v1, the code of the model you are using for vero Beverton-Holt= "BV" , Ricker = "RC", Lotka-Volterra = "LV"
#v2, the code of the model you are using for trcy  
posterior_feasibility<-function(vero_model,trcy_model,vero_growth,trcy_growth,v1,v2){
  
  gi<-.372

  gj<-.258
  #sj<-.033
  #extract the posterior distribution of alphas
  vero_post<-posterior_samples(vero_model)
  trcy_post<-posterior_samples(trcy_model)
  #that should correspond to the same posterior of the _growth models
  r_vero<-vero_growth[,v1]
  r_trcy<-trcy_growth[,v2]
  
  
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(num_posterior){
   
    omega_results       <-c()
    feasibility_results <-c()
    growth_results      <-c()
    theta_results       <-c()

    for( i in 1:nrow(vero_post)){
      
      #we get the corresponding posterior values, vero first, trcy first, gi (vero), gj(trcy)
      alpha  <- get_posterior_alphas(vero_post[i,],trcy_post[i,],gi,gj)

      
       r1     <- r_vero[[i,1]]
       r2     <- r_trcy[[i,1]]
      
       growth <-r2/r1
       #And estimate the feasability domain
       omega       <-Omega(alpha)
       feasibility <-test_feasibility(alpha,c(r1,r2))
       theta       <-theta(alpha,c(r1,r2))
      
      # #we save it 
       omega_results       <-c(omega_results,omega)
       feasibility_results <-c(feasibility_results, feasibility)
       growth_results      <-c(growth_results,growth)
       theta_results       <-c(theta_results,theta)
    }
    
   
    pp<- cbind(omega_results,feasibility_results, growth_results,theta_results)

    pp<-as.data.frame(pp)
    return(pp)
  }else{warning("Posterior distributions are not the same length")}
  
  
}


fixed_feasibility<-function(vero_model,trcy_model,vero_growth,trcy_growth,v1,v2){
  
  gi<-.372
  # si<-.556
  gj<-.258
  #sj<-.033
  params_vero<-fixed_model(vero_model)
  params_trcy<-fixed_model(trcy_model)
  #aii #alphaij
  #aji #alphajj
  alpha<-get_fixed_alphas(params_vero,params_trcy,gi,gj)
  omega<-Omega(alpha)
  #this has to change if you change models!
  r_vero <- vero_growth[,v1]
  r_trcy <- trcy_growth[,v2]
  
  r1<-mean(r_vero[[1]])
  r2<-mean(r_trcy[[1]])
  
  feasibility <-test_feasibility(alpha,c(r1,r2))
  growth      <-r2/r1
  theta       <-theta(alpha,c(r1,r2))
  return(c(omega,feasibility,growth,theta))
}






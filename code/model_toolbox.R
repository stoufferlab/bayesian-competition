#toolbox for example BEV competition
require(brms)
require(tidyverse)
require(ggplot2)
require(rethinking)
source("code/feasibility_toolbox.R")
source("code/extract_params.R")
source("code/growth_rates.R")


get_alphas<-function(vero,trcy){
  vero<-as.list(vero)
  trcy<-as.list(trcy)
  a11<-vero$b_alphaii_Intercept
  a21<-trcy$b_alphaji_Intercept
  a12<-vero$b_alphaij_Intercept
  a22<-trcy$b_alphajj_Intercept
  alpha<-matrix( c(a11,a21,a12,a22) ,ncol=2,nrow=2)
  return(alpha)
}
fixed_coexistence<-function(vero_model,trcy_model,titulo){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  params_vero<-fixed_model(vero_model)
  params_trcy<-fixed_model(trcy_model)
  #aii #alphaij
  #aji #alphajj
  alpha<-matrix(c(params_vero$alphaii_Intercept,params_trcy$alphaji_Intercept,params_vero$alphaij_Intercept,params_trcy$alphajj_Intercept),ncol = 2,nrow=2)
  omega<-Omega(alpha)
  #this has to change if you change models!
  r1<-bev_growth(si,gi,params_vero$lambdai_Intercept)
  r2<-bev_growth(sj,gj,params_trcy$lambdaj_Intercept)
  
  projection_2sp(alpha,r1,r2,titulo)
  feasibility<-test_feasibility(alpha,c(r1,r2))
  
  return(c(omega,feasibility))
  
}

posterior_coexistence<-function(vero_model,trcy_model){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  vero_post<-posterior_samples(vero_model)
  trcy_post<-posterior_samples(trcy_model)
  
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(num_posterior){
  cols<-.01
  omega_results<-c()
  feasibility_results<-c()
  growth_results<-c()
  for( i in 1:nrow(vero_post)){
  
    alpha<-get_alphas(vero_post[i,],trcy_post[i,])
    r1<-bev_growth(si,gi,vero_post[i,]$b_lambdai_Intercept)
    r2<-bev_growth(sj,gj,trcy_post[i,]$b_lambdaj_Intercept)  
    gg<-r2/r1
    results<-projection_posterior(alpha,r1,r2,cols)
    omega_results<-c(omega_results,results[1])
    feasibility_results<-c(feasibility_results, results[2])
    growth_results<-c(growth_results,gg)
  }
 
   vv<-vero_post %>% select(b_lambdai_Intercept, b_alphaii_Intercept, b_alphaij_Intercept )
   tt<-trcy_post  %>% select( b_lambdaj_Intercept, b_alphaji_Intercept, b_alphajj_Intercept)  
   pp<- cbind(omega_results,feasibility_results, growth_results)
   pp<-as.data.frame(pp)
   return(pp)
  }else{warning("Posterior distributions are not the same length")}
  
  
}
posterior_feasibility<-function(vero_model,trcy_model){
  
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  vero_post<-posterior_samples(vero_model)
  trcy_post<-posterior_samples(trcy_model)
  
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(num_posterior){
   
    omega_results       <-c()
    feasibility_results <-c()
    growth_results      <-c()
    theta_results       <-c()
    for( i in 1:nrow(vero_post)){
      
      #we get the corresponding posterior values
      alpha  <-get_alphas(vero_post[i,],trcy_post[i,])
      r1     <-bev_growth(si,gi,vero_post[i,]$b_lambdai_Intercept)
      r2     <-bev_growth(sj,gj,trcy_post[i,]$b_lambdaj_Intercept)  
      growth <-r1/r2
      #And estimate the feasability domain
      omega       <-Omega(alpha)
      feasibility <-test_feasibility(alpha,c(r1,r2))
      theta       <-theta(alpha,c(r1,r2))
     
      #we save it 
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






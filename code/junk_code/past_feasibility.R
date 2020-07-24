#past feasibility functions
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

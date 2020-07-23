fixed_model<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  return(params)
  }


get_posterior_alphas<-function(vero,trcy,gi,gj){
  vero<-as.list(vero)
  trcy<-as.list(trcy)
  a11<-vero$b_alphaii_Intercept
  a21<-trcy$b_alphaji_Intercept
  a12<-vero$b_alphaij_Intercept
  a22<-trcy$b_alphajj_Intercept
  alpha<-matrix( c(a11,a21,a12,a22) ,ncol=2,nrow=2) 
  alpha<-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
  
  return(alpha)
}





get_fixed_alphas<-function(vero,trcy,gi,gj){
  vero<-as.list(vero)
  trcy<-as.list(trcy)
  a11<-vero$alphaii_Intercept
  a21<-trcy$alphaji_Intercept
  a12<-vero$alphaij_Intercept
  a22<-trcy$alphajj_Intercept
  alpha<-matrix( c(a11,a21,a12,a22) ,ncol=2,nrow=2) 
  alpha<- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
  return(alpha)
}
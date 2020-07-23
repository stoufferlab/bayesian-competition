

projection_2sp<- function(alpha,r1,r2,titulo){
  
  #graphics.off()
  #pdf(file=file_name, width = 8, height = 8, useDingbat=FALSE, family = "ArialMT")
  # par(mar = c(3,3,3,3))
  
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  # v3 <- alpha_n[,3]
  # vc <- (v1 + v2 + v3)
  vc <- (v1 + v2 )
  vc <- vc / sqrt(sum(vc^2))
  
  v1 <- v1/sum(v1)
  v2 <- v2/sum(v2)
  #v3 <- v3/sum(v3)
  vc <- vc/sum(vc)
  
  
  
  
  plot(c(0,2),c(0,0),axes=F,col='grey50',type='l',lwd=2,xlim=c(0,1),ylim=c(0,1), main=titulo, xlab="Intrinsic growth rate vero (r1)",ylab="Intrinisc growth rate trcy (r2)")
  lines(c(0,0), c(0,2), col='grey50', lwd=2)
  #####
  
  
  
  vcP <- v1/sqrt(sum(v1^2)) + v2/sqrt(sum(v2^2))
  
  vcP <- vcP / sum(vcP)
  
  
  points( v1[1], v1[2], col='dodgerblue',pch=16,cex=1.5)
  points( v2[1], v2[2], col='dodgerblue',pch=16,cex=1.5)
  points( vcP[1], vcP[2], col="orange3", pch=16,cex=1.5 )
  
  upper<-alpha[2,2]/alpha[1,2]
  lower<-alpha[2,1]/alpha[1,1]
  
  xx<-seq(0,v1[1],.01)
  yy<-seq(0,v2[1],.01)
  lines(xx,xx*lower,col='mediumseagreen',lty=2,lwd=2)
  lines(yy,yy*upper,col='mediumseagreen',lty=2,lwd=2)
  
  growth<-r2/r1
  zz<-seq(0,vcP[1],0.01)
  lines(zz,zz*growth, col="brown",lwd=2)
  
  centroid<-r_centroid(alpha)
  points(centroid[1,],centroid[2,],col="orange1",pch=18)
  
}

projection_posterior<-function(alpha,r1,r2,col_alpha){
  omega<-Omega(alpha)
  feasibility<-test_feasibility(alpha,c(r1,r2))
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  vc <- (v1 + v2 )
  vc <- vc / sqrt(sum(vc^2))
  v1 <- v1/sum(v1)
  v2 <- v2/sum(v2)
  vc <- vc/sum(vc)
  vcP <- v1/sqrt(sum(v1^2)) + v2/sqrt(sum(v2^2))
  vcP <- vcP / sum(vcP)
  
  cblue<-rethinking::col.alpha("dodgerblue",alpha = col_alpha)
  corange<-rethinking::col.alpha("orange3",alpha = col_alpha)
  cgreen<-rethinking::col.alpha("mediumseagreen",alpha = col_alpha)
  cbrown<-rethinking::col.alpha("brown",alpha = col_alpha)
  
  points( v1[1], v1[2], col=cblue,pch=16,cex=1.5)
  points( v2[1], v2[2], col=cblue,pch=16,cex=1.5)
  points( vcP[1], vcP[2], col=corange, pch=16,cex=1.5 )
  
  upper<-alpha[2,2]/alpha[1,2]
  lower<-alpha[2,1]/alpha[1,1]
  
  xx<-seq(0,v1[1],.01)
  yy<-seq(0,v2[1],.01)
  lines(xx,xx*lower,col=cgreen,lty=2,lwd=1)
  lines(yy,yy*upper,col=cgreen,lty=2,lwd=1)
  
  growth<-r2/r1
  zz<-seq(0,vcP[1],0.01)
  lines(zz,zz*growth, col=cbrown,lwd=2)
  
  return(c(omega,feasibility))
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
  theta      <- theta(alpha,c(r1,r2))
  growth     <- r2/r1
  return(c(omega,feasibility,theta,growth))
  
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

#######fixed two species

#toolbox for example BEV competition
require(brms)
require(tidyverse)
require(ggplot2)
require(rethinking)
source("code/feasibility_toolbox.R")
source("code/model_toolbox.R")
source("code/growth_rates.R")

fixed_model<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  return(params)
}

alpha_matrix_fixed <- function(vero_pars, trcy_pars, gi,gj, env){
  # vero_pars <- as.list(vero_pars)
  # trcy_pars <- as.list(trcy_pars)
  
  if(env){
    alpha11 <- vero_pars$alphaii_Intercept + vero_pars$alphaii_env
    alpha21 <- trcy_pars$alphaij_Intercept + trcy_pars$alphaij_env
    alpha12 <- vero_pars$alphaij_Intercept + vero_pars$alphaij_env
    alpha22 <- trcy_pars$alphaii_Intercept + trcy_pars$alphaii_env
    alpha   <-matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2)
    alpha   <-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }else{
    alpha11 <- vero_pars$alphaii_Intercept
    alpha21 <- trcy_pars$alphaij_Intercept
    alpha12 <- vero_pars$alphaij_Intercept
    alpha22 <- trcy_pars$alphaii_Intercept
    alpha   <-matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2) 
    alpha   <-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }
  return(alpha)  
}


alpha_matrix <- function(vero_row, trcy_row, gi,gj, env){
  # vero_row <- as.list(vero_row)
  # trcy_row <- as.list(trcy_row)
  
  if(env){
    alpha11 <- vero_row$b_alphaii_Intercept + vero_row$b_alphaii_env
    alpha21 <- trcy_row$b_alphaij_Intercept + trcy_row$b_alphaij_env
    alpha12 <- vero_row$b_alphaij_Intercept + vero_row$b_alphaij_env
    alpha22 <- trcy_row$b_alphaii_Intercept + trcy_row$b_alphaii_env
    alpha   <-matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2)
    alpha   <-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }else{
    alpha11 <- vero_row$b_alphaii_Intercept
    alpha21 <- trcy_row$b_alphaij_Intercept
    alpha12 <- vero_row$b_alphaij_Intercept
    alpha22 <- trcy_row$b_alphaii_Intercept
    alpha   <-matrix( c(alpha11,alpha21,alpha12,alpha22) ,ncol=2,nrow=2) 
    alpha   <-sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
    
  }
  return(alpha)  
}







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
  
  
  
  
  plot(c(0,2),c(0,0),axes=F,col='grey50',type='l',lwd=2,xlim=c(0,1),ylim=c(0,1), main=titulo, xlab="Intrinsic growth rate vero",ylab="Intrinisc growth rate trcy")
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


feasibility_params<-function(vero_model,vero_fun, trcy_model, trcy_fun){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  vero_pars<-fixed_model(BEV_vero)
  trcy_pars<-fixed_model(BEV_trcy)
  
  alpha<- alpha_matrix_fixed(vero_pars = vero_pars, trcy_pars = trcy_pars, gi=gi, gj=gj, env=0)
  
  vero_growth<-vero_fun(si,gi,vero_pars$lambda_Intercept)
  trcy_growth<-trcy_fun(sj,gj,trcy_pars$lambda_Intercept)
  
  params <-list("alpha"=alpha,"growth_vero"=vero_growth,"growth_trcy"=trcy_growth)
  
  return(params)
  
}


test<-feasibility_params(vero_model = BEV_vero, vero_fun = bev_growth, trcy_model = BEV_trcy, trcy_fun = bev_growth)

omega<- Omega(test$alpha)
theta_r<-theta(test$alpha, c(test$growth_vero,test$growth_trcy))
test_feasibility(test$alpha, c(test$growth_vero,test$growth_trcy))

projection_2sp(test$alpha,test$growth_vero,test$growth_trcy,"Beverton-Holt")


vero_post<-posterior_parameters(model=BEV_vero,fun = bev_growth,si,gi,0)
trcy_post<-posterior_parameters(model=BEV_trcy,fun = bev_growth,sj,gj,0)


vero_mean<-as.list(apply(vero_post,2,mean))
trcy_mean<-as.list(apply(trcy_post,2,mean))








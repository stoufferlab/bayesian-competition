## theoretical figure

source("code/growth_rates.R")
source("code/feasibility_toolbox.R")
source("code/model_toolbox.R")

vero_bh_multispecies_poisson <- readRDS("~/bayesian-competition/model_objects/vero_bh_multispecies_poisson.rds")
trcy_bh_multispecies_poisson <- readRDS("~/bayesian-competition/model_objects/trcy_bh_multispecies_poisson.rds")

vero_model<- vero_bh_multispecies_poisson
trcy_model<- trcy_bh_multispecies_poisson

gi<-.372
si<-.556
gj<-.258
sj<-.033

fixed_model<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  return(params)
}

get_fixed_alphas<-function(vero_model,gi,trcy_model,gj){
  vero<-fixed_model(vero_model)
  trcy<-fixed_model(trcy_model)
  
  a11<-vero$alphaii_Intercept
  a21<-trcy$alphaij_Intercept
  a12<-vero$alphaij_Intercept
  a22<-trcy$alphaii_Intercept
  
  alpha<-matrix( c(a11,a21,a12,a22) ,ncol=2,nrow=2) 
  alpha<- sweep (alpha,MARGIN=2,STAT=c(gi,gj),FUN="*")
  return(alpha)
}

alpha<- get_fixed_alphas(vero_model = vero_bh_multispecies_poisson, gi = gi, trcy_model = trcy_bh_multispecies_poisson, gj = gj)

get_fixed_growth<- function(vero_model,vero_function, si,gi, trcy_model, trcy_function, sj,gj){
  vero<-fixed_model(vero_model)
  trcy<-fixed_model(trcy_model)
  
  
  r_vero<- vero_function(s=si,g=gi, lambda= exp(vero$lambda_Intercept))
  r_trcy<- trcy_function(s=sj,g=gj, lambda= exp(trcy$lambda_Intercept)) 
  return(c(r_vero,r_trcy))
}

rr<- get_fixed_growth(vero_model = vero_bh_multispecies_poisson,vero_function = bh_growth,si = si,gi = gi,
                      trcy_model = trcy_bh_multispecies_poisson,trcy_function = bh_growth,sj = sj,gj = gj)



projection_2sp<- function(alpha,r1,r2, title, sp1=FALSE, sp2 =FALSE){
  
  if(sp1){ xlb= "Intrinsic growth rate sp. 1"}else{xlb= ""}
  
  if(sp2){ ylb= "Intrinsic growth rate sp. 2"}else{ylb= ""}
  
  
  
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[, 1]
  v2 <- alpha_n[, 2]
  vc <- (v1 + v2)
  vc <- vc / sqrt(sum(vc ^ 2))
  
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  #v3 <- v3/sum(v3)
  vc <- vc / sum(vc)
  
  
  
  
  plot(
    c(0, 2),
    c(0, 0),
    axes = F,
    col = 'grey50',
    type = 'l',
    lwd = 2,
    xlim = c(0, .8),
    ylim = c(0, .8),
    xlab = xlb,
    ylab = ylb,
    cex.lab = 1.5,
    cex.main = 2,
    main = title,
    adj  = 0
  )
  lines(c(0, 0), c(0, 2), col = 'grey50', lwd = 2)
  #####
  
  
  
  vcP <- v1 / sqrt(sum(v1 ^ 2)) + v2 / sqrt(sum(v2 ^ 2))
  
  vcP <- vcP / sum(vcP)
  
  
  points( v1[1], v1[2], col='mediumseagreen',pch=16,cex=1)
  points( v2[1], v2[2], col='mediumseagreen',pch=16,cex=1)
  points( vcP[1], vcP[2], col="orange1", pch=16,cex=1)
  
  
  
  
  upper <- alpha[2, 2] / alpha[1, 2]
  lower <- alpha[2, 1] / alpha[1, 1]
  
  xx <- seq(0, v1[1], .01)
  yy <- seq(0, v2[1], .01)
  lines(xx,
        xx * lower,
        col = 'mediumseagreen',
        lty = 1,
        lwd = 2)
  lines(yy,
        yy * upper,
        col = 'mediumseagreen',
        lty = 1,
        lwd = 2)
  
  
  
  growth <- r2 / r1
  zz <- seq(0, vcP[1]-.01, 0.01)
  # lines(zz,
  #       zz * growth,
  #       col = "tomato3",
  #       lwd = 1.5,
  #       lty = 2)
  arrows(0,0, max(zz), max(zz*growth), col="tomato3", lwd=2)
  
  print(max(zz*growth))
  centroid <- r_centroid(alpha)
  #points(centroid[1, ], centroid[2, ], col = "orange1", pch = 18)
  
  slope_centroid <- (centroid[2, ] - vcP[2]) / (centroid[1, ] - vcP[1])
  
  lines(zz, slope_centroid * zz, lwd = 2, col = "orange1")
  
 
  
}




move_alpha<-function(alpha, x){
  alpha[1,1]<- alpha[1,1] -x
  alpha[1,2]<- alpha[1,2] + x
  
  alpha[2,1] <- alpha[2,1] + x
  alpha[2,2] <- alpha[2,2] - x
  
  return(alpha)
}



moved_projection<-function(alpha, r1, r2,x=0.005, title, sp1 = FALSE,sp2 = FALSE){
   nn <-rnorm(n = 500,mean = 0,sd = .0015)
 #  nn <- seq(0,x, 0.0001)
   rr <- rnorm(n = length(nn),mean = 0, sd = 0.7)
   
   projection_2sp(alpha = alpha, r1 = r1, r2 = r2, title = title, sp1=sp1, sp2=sp2)
   
   for(i in 1:length(nn)){
     col_g <-rethinking::col.alpha(acol = 'mediumseagreen' , alpha = 0.02)
     col_o <-rethinking::col.alpha(acol = 'orange1' , alpha = 0.02)
     col_b <-rethinking::col.alpha(acol = 'tomato3' , alpha = 0.01)
     
     
     alpha_new <- move_alpha(alpha = alpha, x=nn[i])
     
     D <- diag(1 / sqrt(diag(t(alpha_new) %*% alpha_new)))
     alpha_n <- alpha_new %*% D
     
     v1 <- alpha_n[, 1]
     v2 <- alpha_n[, 2]
     vc <- (v1 + v2)
     vc <- vc / sqrt(sum(vc ^ 2))
     v1 <- v1 / sum(v1)
     v2 <- v2 / sum(v2)
     vc <- vc / sum(vc)
     
     vcP <- v1 / sqrt(sum(v1 ^ 2)) + v2 / sqrt(sum(v2 ^ 2))
     vcP <- vcP / sum(vcP)
     
     points( v1[1], v1[2], col= col_g  ,pch=16,cex=1)
     points( v2[1], v2[2], col= col_g, pch=16,cex=1)
     points( vcP[1], vcP[2], col=col_o, pch=16,cex=1)
     
       growth <- (r2 + rr[i]) / (r1  + ( rr[i]))
      
     max_h<-0.7394366
     zh <- max_h/growth
     
      arrows(0,0, zh, max_h, col=col_b, lwd=2)
      
    
     
   }
}

no_coexist<-matrix( c( 0.02,0.023,0.025,0.015), 2,2)



pdf(file="results/theory_figure.pdf", width = 7.7, height = 7.7/3 )
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE),
       widths=c(1,1,1), heights=c(1,1,1))

x<-1.5
par(mgp=c(.5,.5,0), mar=c(x,x,x,x))
projection_2sp(alpha = alpha, r1 = rr[1], r2 = 9, title = "A",sp1 = FALSE,sp2 = TRUE)
moved_projection(alpha = alpha, r1 = 4.26, r2=9, x= 0.005, title ="B",sp1 = TRUE,sp2 = FALSE)
projection_2sp(alpha = no_coexist, r1 = 3, r2 = 4, title = "C", sp1 = FALSE,sp2 = FALSE)

dev.off()

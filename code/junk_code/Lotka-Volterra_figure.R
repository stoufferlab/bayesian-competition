gi<-.372
si<-.556
gj<-.258
sj<-.033

Ni_isocline<-function(si,gi,lambdai,alphaii,alphaij,Nj){
seeds<--si*(1-gi)
competition<- lambdai* (1-(alphaij*Nj))
denominator<- alphaii*lambdai
Ni<-(seeds+competition+1)/denominator
Ni
}
Nj_isocline<-function(sj,gj,lambdaj,alphajj,alphaji,Ni){s
  seeds<--sj*(1-gj)
  competition<- lambdaj* (1-(alphaji*Ni))
  denominator<- alphajj*lambdaj
  Nj<-(seeds+competition+1)/denominator
  Nj
}
isoclines<-function(max_seeds,si,gi,lambdai,alphaii,alphaij,sj,gj,lambdaj,alphajj,alphaji,col1,col2) {
  N_x<-seq(0,500)
  Ni<- Ni_isocline(si,gi,lambdai,alphaii,alphaij,N_x)
  Nj<-Nj_isocline(sj,gj,lambdaj,alphajj,alphaji,N_x)
  lines(Ni,N_x,lwd=2,col=col1)
  lines(N_x,N_j,lwd=2, col=col2)
}


colour_generator<-function(red,green,blue ,alpha){
  r<-red/255
  g<-green/255
  b<-blue/255
  col<-rgb(r,g,b,alpha)
}


parnames1<-c("b_lambdai_Intercept","b_alphaii_Intercept","b_alphaij_Intercept")
parnames2<-c("b_lambdaj_Intercept","b_alphajj_Intercept", "b_alphaji_Intercept")
posterior_points_alphas<-function(model1,parnames1,model2,parnames2,n_samples){
  sub<-seq(1,n_samples,1)
  post1<-posterior_samples(model1,parnames1,exact_match = T, subset=sub)
  lambdai<-post1[,1]
  alphaii<-post1[,2]
  alphaij<-post1[,3]
  post2<-posterior_samples(model2,parnames2,exact_match = T, subset=sub)
  lambdaj<-post2[,1]
  alphajj<-post2[,2]
  alphaji<-post2[,3]
  points_refill<-data.frame("lambdai"=lambdai, "alphaii"=alphaii,"alphaij"=alphaij,"lambdaj"=lambdaj , "alphajj"=alphajj, "alphaji"=alphaji )
  points_refill
}
post<-posterior_points_alphas(LV_vero,parnames1,LV_trcy,parnames2,1000)

zngi_posterior<-function(post,si,gi,sj,gj){
 nn<-nrow(post)
  

  
  col1<-colour_generator(41,24,164,1)
  col1.5<-colour_generator(41,24,164,.01)
  col2<-colour_generator(147,164,24,1)
  col2.5<-colour_generator(147,164,24,.01)
  
  plot(NULL, xlim=c(0,500), ylim=c(0,500), xlab="Ni", ylab="Nj")
  mean_post<-apply(post,2,mean)
  isoclines(500,si,gi,mean_post[1],mean_post[2], mean_post[3],sj,gj, mean_post[4], mean_post[5],mean_post[6],col1,col2)
  
  for(i in 1:nn){
    isoclines(500,si,gi,post[i,1],post[i,2],post[i,3],sj,gj,post[i,4],post[i,5],post[i,6],col1.5,col2.5)
  }
}
zngi_posterior(post,si,gi,sj,gj)


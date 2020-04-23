gi<-.372
si<-.556
gj<-.258
sj<-.033
post_exp<-posterior_alphas(exp_vero,parnames1,exp_trcy,parnames2,8000)



Ni_invader_e<-posterior_Ni_invader_exp(post_exp,gi,si,gj,sj)
Nj_invader_e<-posterior_Nj_invader_exp(post_exp,gi,si,gj,sj)

coexistence_exp<-function(Ni_invasion,Nj_invasion){
  nn<-length(Ni_invasion)
  outcome<-c()
  for(i in 1:nn){
    outcome[i]<-ifelse(Ni_invasion[i]>0 && Nj_invasion[i]>0,1,0)
  }
  predictions<-data.frame("Ni"=Ni_invasion, "Nj"=Nj_invasion, "Coexistence"=outcome)
  predictions
}
mutual_invasibility<-coexistence_exp(Ni_invader_e, Nj_invader_e)
mutual_invasibility$Coexistence

colour_generator<-function(red,green,blue ,alpha){
  r<-red/255
  g<-green/255
  b<-blue/255
  col<-rgb(r,g,b,alpha)
}

coexistence_graph<-function(mutual_invasibility){
  plot(NULL,xlab="Logarithmic invastion growth rate of vero", ylab="Logarithmic invasion growth rate of tryc", xlim=c(-1,3),ylim=c(-1,3), main="Exponential model")
  
  lines(c(0,0),c(-10,10))
  lines(c(-10,10),c(0,0))
  nn<-nrow(mutual_invasibility)
  col1<-colour_generator(24,164,151,.3)
  col2<-colour_generator(164,24,37,.3)
  pp<-apply(post_exp,2,mean, na.rm=T)
  #lambdai,alphaij,lambdaj,alphajj,gi,si,gj,sj
  mean_vero<-Ni_invader_exp(pp[1],pp[3],pp[4],pp[5],gi,si,gj,sj)
  #lambdai,alphaji,lambdaj,alphaii,gi,si,gj,sj
  mean_trcy<-Nj_invader_exp(pp[1],pp[6],pp[4],pp[2], gi,si,gj,sj)
  points(mean_vero,mean_trcy, pch=17, col=1, cex=1.5)
  nn<-nrow(mutual_invasibility)
  col1<-colour_generator(24,164,151,.3)
  col2<-colour_generator(164,24,37,.3)
  for(i in 1:nn){
    if(mutual_invasibility[i,3]==1){
      points(mutual_invasibility[i,1], mutual_invasibility[i,2],pch=20,col=col1)
    }else{
      points(mutual_invasibility[i,1], mutual_invasibility[i,2],pch=20,col=col2)
    }
  }
  points(mean_vero,mean_trcy, pch=17, col="goldenrod1", cex=1)
}

coexistence_graph(mutual_invasibility)



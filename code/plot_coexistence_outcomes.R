##plot invasion growth rates
gi<-.372
si<-.556
gj<-.258
sj<-.033
post<-posterior_alphas(bev_vero,parnames1,bev_trcy,parnames2,8000)

Ni_invasion<-posterior_Ni_invader(post,gi,si,gj,sj)
Nj_invasion<-posterior_Nj_invader(post,gi,si,gj,sj)

coexistence<-function(Ni_invasion,Nj_invasion){
  nn<-length(Ni_invasion)
  outcome<-c()
  for(i in 1:nn){
    outcome[i]<-ifelse(Ni_invasion[i]>1 && Nj_invasion[i]>1,1,0)
  }
  predictions<-data.frame("Ni"=Ni_invasion, "Nj"=Nj_invasion, "Coexistence"=outcome)
  predictions
}
mutual_invasibility<-coexistence(Ni_invasion,Nj_invasion)
mutual_invasibility$Coexistence

colour_generator<-function(red,green,blue ,alpha){
  r<-red/255
  g<-green/255
  b<-blue/255
  col<-rgb(r,g,b,alpha)
}

coexistence_graph<-function(mutual_invasibility,gi,si,gj,sj){
  plot(NULL,xlab="Invastion growth rate of vero", ylab="Invasion growth rate of tryc", xlim=c(-1,3),ylim=c(-1,3), main="Beverton Holt model")
  lines(c(1,1),c(-10,10))
  lines(c(-10,10),c(1,1))
  nn<-nrow(mutual_invasibility)
  col1<-colour_generator(24,164,151,.3)
  col2<-colour_generator(164,24,37,.3)
  pp<-apply(post,2,mean, na.rm=T)
  #lambdai,alphaij,lambdaj,alphajj,gi,si,gj,sj
   mean_vero<-Ni_invader(pp[1],pp[3],pp[4],pp[5],gi,si,gj,sj)
   #lambdai,alphaji,lambdaj,alphaii,gi,si,gj,sj
   mean_trcy<-Nj_invader(pp[1],pp[6],pp[4],pp[2], gi,si,gj,sj)
   points(mean_vero,mean_trcy, pch=17, col=1, cex=1.5)
  for(i in 1:nn){
    if(mutual_invasibility[i,3]==1){
      points(mutual_invasibility[i,1], mutual_invasibility[i,2],pch=20,col=col1)
    }else{
      points(mutual_invasibility[i,1], mutual_invasibility[i,2],pch=20,col=col2)
    }
  }
   points(mean_vero,mean_trcy, pch=17, col="goldenrod1", cex=1)
   

}

coexistence_graph(mutual_invasibility,gi,si,gj,sj)









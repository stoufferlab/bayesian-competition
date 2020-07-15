##calculate growth rate when rare
library(brms)
library(dplyr)

parnames1<-c("b_lambdai_Intercept","b_alphaii_Intercept","b_alphaij_Intercept")
parnames2<-c("b_lambdaj_Intercept","b_alphajj_Intercept", "b_alphaji_Intercept")
posterior_alphas<-function(model1,parnames1,model2,parnames2, n_samples){
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



Ni_invader<-function(lambdai,alphaij,lambdaj,alphajj,gi,si,gj,sj){

  Nj_equilibrium<-(1/ (alphajj))* ( ( (lambdaj) / (1 - ((1-gj)*sj) )  )  -1)
  fi<-(lambdai)/(1 + (alphaij*Nj_equilibrium))
  Ni_growth<-((1-gi)*si) + (fi)
  Ni_growth
}

Nj_invader<-function(lambdai,alphaji,lambdaj,alphaii,gi,si,gj,sj){
  Ni_equilibrium<-(1/ (alphaii))* ( ( (lambdai) / (1 - ((1-gi)*si) )  )  -1)
  fj<-(lambdaj)/(1 + (alphaji*Ni_equilibrium))
  Nj_growth<-((1-gj)*sj) + (fj)
  Nj_growth
}


posterior_Ni_invader<-function(post,gi,si,gj,sj){
  pars<-post %>% select(lambdai,alphaij,lambdaj,alphajj)
  nn<-nrow(post)
  Ni<-c()
  for(i in 1:nn){
    Ni[i]<-Ni_invader(pars[i,1], pars[i,2],pars[i,3],pars[i,4],gi,si,gj,sj)
  }
  Ni
}
posterior_Nj_invader<-function(post,gi,si,gj,sj){
  pars<-post %>% select(lambdai,alphaji,lambdaj,alphaii)
  nn<-nrow(post)
  Nj<-c()
  for(i in 1:nn){
    Nj[i]<-Nj_invader(pars[i,1], pars[i,2],pars[i,3],pars[i,4],gi,si,gj,sj)
  }
  Nj
}

#Ni_invasion_bev<-posterior_Ni_invader(post,.372,.556)
#Nj_invasion_bev<-posterior_Nj_invader(post,.258,.033)
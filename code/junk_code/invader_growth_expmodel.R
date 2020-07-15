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



Ni_invader_exp<-function(lambdai,alphaij,lambdaj,alphajj,gi,si,gj,sj){
  Nj_equilibrium<- log( ((-sj*(1-gj)) + 1 )/lambdaj  ) / ((-alphajj)/gj)
  fi<-lambdai*exp((-alphaij/gj)*Nj_equilibrium)
  Ni_growth<-((1-gi)*si) + (fi)
  Ni_growth<-log(Ni_growth)
  Ni_growth
}



Nj_invader_exp<-function(lambdai,alphaji,lambdaj,alphaii,gi,si,gj,sj){
  Ni_equilibrium<- log( ((-si*(1-gi)) + 1 )/lambdai  ) / ((-alphaii)/gi)
  fj<-lambdaj*exp((-alphaji/gi)*Ni_equilibrium)
  Nj_growth<-((1-gj)*sj) + (fj)
  Nj_growth<-log(Nj_growth)
  Nj_growth
}




posterior_Ni_invader_exp<-function(post,gi,si,gj,sj){
  pars<-post %>% select(lambdai,alphaij,lambdaj,alphajj)
  nn<-nrow(post)
  Ni<-c()
  for(i in 1:nn){
    Ni[i]<-Ni_invader_exp(pars[i,1], pars[i,2],pars[i,3],pars[i,4],gi,si,gj,sj)

  }
 Ni
}




posterior_Nj_invader_exp<-function(post,gi,si,gj,sj){
  pars<-post %>% select(lambdai,alphaji,lambdaj,alphaii)
  nn<-nrow(post)
  Ni<-c()
  for(i in 1:nn){
    Ni[i]<-Ni_invader_exp(pars[i,1], pars[i,2],pars[i,3],pars[i,4],gi,si,gj,sj)
  }
  Ni
}



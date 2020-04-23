

cc<-function(m1,m2, samples){
  
  
  m1<-add_criterion(m1, "waic",nsamples=samples)
  m2<-add_criterion(m2, "waic",nsamples=samples)
  
  w<-loo_compare(m1,m2,criterion="waic")
  print(w)
  m<-model_weights(m1,m2, weights="waic")
  print(m)
  
  
  
}

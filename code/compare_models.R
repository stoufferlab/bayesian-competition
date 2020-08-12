

cc<-function(m1,m2,m3,m4){
  
  
  m1<-add_criterion(m1, "waic")
  m2<-add_criterion(m2, "waic")
 m3<-add_criterion(m3, "waic")
 m4<-add_criterion(m4, "waic")
  w<-loo_compare(m1,m2,m3,m4,criterion="waic")
  print(w)
  m<-model_weights(m1,m2,m3,m4,weights="waic")
  print(m)
  
  
  
}

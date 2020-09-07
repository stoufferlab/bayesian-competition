

cc<-function(m1,m2,m3,m4,m5){
  
  
  m1<-add_criterion(m1, "waic")
  m2<-add_criterion(m2, "waic")
m3<-add_criterion(m3, "waic")
 m4<-add_criterion(m4, "waic")
 m5<-add_criterion(m4, "waic")
  w<-loo_compare(m1,m2,m3,m4,m5,criterion="waic")
  print(w)
  m<-model_weights(m1,m2,m3,m4,m5,weights="waic")
  print(m)
  
  
  
}


model_compare<-function(...){
   models <- list(...)
   waic_value <- c()
   looic_value <-c()
   
   for(i in 1:length(models)){
     w <- waic(models[[i]])
     wa <- w$waic
     waic_value[i] <-wa
     
     l <- loo(models[[i]])
     lo <- l$looic
     looic_value[i] <-lo
   }
  
  waic_weights <- model_weights(..., weights = "waic")
  loo_weights <- model_weights(..., weights = "loo")
  out <- cbind(waic_value, waic_weights, looic_value, loo_weights)
  out <- out[order(waic_weights),]
  return(out)
  
}

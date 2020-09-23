#Function to compare models, you give it as many brms objects as you want!




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

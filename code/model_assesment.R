library(ggplot2)
library(ggridges)
library(dplyr)
library(bayesplot)
library(tidybayes)
        
        
convergence<-function(model){
  post      <-posterior_samples(model, add_chain = TRUE)
  trace     <-mcmc_trace(post) 
  densities <-mcmc_dens(post)
  rhat_vals <- rhat(model)
 # mcmc_rhat_data(rhat_vals) 
  rhat_plot       <-mcmc_rhat(rhat_vals) + theme_bw()
#  autocorrelation <-mcmc_acf(post)
  
  
  print(trace)
  print(densities)
  print(rhat_plot)
 # print(autocorrelation)
  
  print(pp_check(model))
  print(pp_check(model, type="rootogram"))
}      
deviance_residuals<-function (y,mu,wt) {
 
  r <- mu * wt
  p <- which(y > 0)
  r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
  2 * r
}

prediction_quantiles<- function(model, data){
  
   
   predictions_data   <- data %>% add_predicted_draws(model)  %>% mutate( residuals = totalseeds - .prediction) %>% mutate( stand_residuals = residuals/.prediction)
   
   predictions_data$totalcompetition <- predictions_data$conspecifics + predictions_data$heterospecifics
   
   deviance_data      <- deviance_residuals( y=predictions_data$totalseeds, mu=predictions_data$.prediction, wt=1) 
   
   
   all_residuals  <- data.frame("totalseeds"=predictions_data$totalseeds,"predictions"= predictions_data$.prediction, "residuals"= predictions_data$residuals,
                                "deviance_residuals" = deviance_data, "stan_res" = predictions_data$stand_residuals)
   
   # p<- ggplot(preds, aes(y=totalseeds, x=total_competition, shape = as.factor(env))) +
   #   geom_point(col= rethinking::col.alpha("dodgerblue",.7)) +
   #   theme_bw() +
   #   stat_pointinterval(aes(x=total_competition,y=.prediction),  col= rethinking::col.alpha("grey50", alpha = 0.5) )
   # 
   # p2 <- ggplot(preds)+
   #   stat_pointinterval(aes(x=totalseeds, y=.prediction), col="dodgerblue2") + 
   #   theme_bw() + 
   #   ylim(0,40) +
   #   geom_abline(slope = 1, intercept = 0, col= "grey50") 
   # 
    p3 <-  ggplot(predictions_data) + 
       stat_pointinterval(aes(x=totalseeds, y=.prediction), col="grey50") +
       # stat_pointinterval(aes(x=totalseeds, y= stan_res), col="dodgerblue2") +
       # stat_pointinterval(aes(x=totalseeds, y= deviance_residuals), col="firebrick") +
        theme_bw() +
        ylim(-10,30)
  

   #print(p)
# #print(p2)
print(p3)

   
}


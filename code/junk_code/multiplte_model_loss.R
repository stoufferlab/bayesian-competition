source("code/read_models.R")
source("code/per_capita_loss.R")
source("code/model_toolbox.R")
source("code/growth_rates.R")

#for the fixed effects of models
fixed_params<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  params
}

multiple_loss <-function(models, model_names, loss_fucntions, neighbors) {
    
  N <- seq(1, neighbors, 1)
  all_loss<- data.frame()
    
    
    for (i in 1:length(models)) {
      
      params <-fixed_params(models[[i]])

      
      if(model_names[i]=="Hassell"){
        loss <- loss_fucntions[[i]]
        conspecific_response <-loss(alpha = params$alphaii_Intercept, N = N, beta = params$beta_Intercept)
        heterospecific_response <-loss(alpha= params$alphaij_Intercept, N=N, beta=params$beta_Intercept)
        response<- as.data.frame(cbind(N, conspecific_response, heterospecific_response))
        response$model <- model_names[i]
        
      }else{
        loss <- loss_fucntions[[i]]
        conspecific_response <-loss(alpha = params$alphaii_Intercept, N = N)
        heterospecific_response <-loss(alpha= params$alphaij_Intercept, N=N)
        response<- as.data.frame(cbind(N, conspecific_response, heterospecific_response))
        response$model <- model_names[i]
        
        
       
      }
      
        
    
        all_loss<- rbind(all_loss, response)
    }
    
  return(all_loss)
    
}



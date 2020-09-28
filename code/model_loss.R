



#for the fixed effects of models
fixed_params<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  params
}


# To calculate the per capita loss of many models
multiple_loss <-function(models, model_names, loss_fucntions, neighbors) {
    
  N <- seq(.1, neighbors, .1)
  all_loss<- data.frame()
  
    for (i in 1:length(models)) {
      
      params <-fixed_params(models[[i]])

      if(model_names[i]=="Hassell"){
        loss <- loss_fucntions[[i]]
      
          conspecific_response <- loss(
          alpha = params$alphaii_Intercept,
          N = N,
          beta = exp(params$beta_Intercept)
        )
        
        conspecific_response_env <-loss(
            alpha = (params$alphaii_Intercept + params$alphaii_env),
            N = N,
            beta = exp(params$beta_Intercept + params$beta_env)
          )
        
        
        heterospecific_response <- loss(
          alpha = params$alphaij_Intercept,
          N = N,
          beta = exp(params$beta_Intercept)
        )
       
        heterospecific_response_env <- loss(
          alpha = (params$alphaij_Intercept + params$alphaij_env),
          N = N,
          beta = exp(params$beta_Intercept + params$beta_env)
        )
        
      
        response <-
          as.data.frame(
            cbind(
              N,
              conspecific_response,
              conspecific_response_env,
              heterospecific_response,
              heterospecific_response_env
            )
          )
        response$model <- model_names[i]
        
      }else{
      
        loss <- loss_fucntions[[i]]
        
        conspecific_response <- loss(
          alpha = params$alphaii_Intercept,
          N = N)
        
        conspecific_response_env <-loss(
          alpha = (params$alphaii_Intercept + params$alphaii_env),
          N = N)
        
        
        heterospecific_response <- loss(
          alpha = params$alphaij_Intercept,
          N = N)
        
        heterospecific_response_env <- loss(
          alpha = (params$alphaij_Intercept + params$alphaij_env),
          N = N)
        
        
        response <-
          as.data.frame(
            cbind(
              N,
              conspecific_response,
              conspecific_response_env,
              heterospecific_response,
              heterospecific_response_env
            )
          )
        response$model <- model_names[i]
        
       
      }
      
        
    
        all_loss<- rbind(all_loss, response)
    }
    
  return(all_loss)
    
}



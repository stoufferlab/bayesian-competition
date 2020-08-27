library(brms)
source("code/models/set_priors.R")
source("code/clean_data.R")

#the model to fit
beverton_holt<-bf(
  totalseeds~ ( exp(lambda) /(1 + (alphaii*conspecifics) + (alphaij*heterospecifics)) ) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)


model <- bf(
  totalseeds ~ exp(lambda -  log(1 + (alphaii*conspecifics) + (alphaij*heterospecifics)) ),
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl = TRUE
)

linear <- bf(
  totalseeds ~ 1 + conspecifics + heterospecifics)





#vero as focal
BEV_vero_poisson<-brm( formula = beverton_holt,
               prior   = prior,
               data    = vero_focal,
               family  = poisson(),
               iter    = 20000,
               warmup  = 16000, 
               cores   = 4,
               chains  = 4,
               inits   = 0 ,
               control =list(adapt_delta = .99, max_treedepth=13)
)



BEV_vero_model<-brm( formula = model,
                       prior   = prior,
                       data    = vero_focal,
                       family  = poisson(),
                       iter    = 20000,
                       warmup  = 16000, 
                       cores   = 4,
                       chains  = 4,
                       inits   = 0 ,
                       control =list(adapt_delta = .99, max_treedepth=13)
)


BEV_vero_linear<-brm( formula = linear,
              
                     data    = vero_focal,
                     family  = poisson(),
                     iter    = 20000,
                     warmup  = 16000, 
                     cores   = 4,
                     chains  = 4,
                     inits   = 0 ,
                     control =list(adapt_delta = .99, max_treedepth=13)
)

saveRDS(BEV_vero,file="BEV_vero.RDS")
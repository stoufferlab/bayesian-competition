#Ricker models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")
source("code/clean_data.R")

#the ricker model
ricker<- bf(
  totalseeds~ (exp(lambda) *exp( -(alphaii*conspecifics)  - (alphaij*heterospecifics)))  ,
       lambda  ~ 1 + env,
       alphaii ~ 1 + env,
       alphaij ~ 1 + env,
       nl=T)
#with vero as focal
RC_vero<-brm(formula = ricker,
              prior  = prior,
              data   = vero_focal,
             family  = poisson(link="identity"),
              iter   = 20000,
             warmup  = 16000, 
             inits   = 0, 
             cores   = 4, 
             chains  = 4,
             control = list(adapt_delta = .99, max_treedepth=13)
)

saveRDS(RC_vero,file="RC_vero.RDS")


#with trcy as focal
RC_trcy<-brm( formula = ricker,
              prior   = prior,
              data    = trcy_focal,
              family  = poisson(link="identity"),
              iter    = 20000,
              warmup  = 16000,
              cores   = 4 ,
              chains  = 4,
              inits   = 0,
              control = list(adapt_delta = .99, max_treedepth=13)
)

saveRDS(RC_trcy,file="RC_trcy.RDS")

#LV models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")
source("code/clean_data.R")

lotka_volterra<-bf(
  totalseeds~ ( exp(lambda) * (1 - (alphaii*conspecifics) - (alphaij*heterospecifics)) ) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)




LV_vero<-brm(formula = lotka_volterra,
             prior   = prior,
             data    = vero_focal,
             family  = poisson("identity"),
             iter    = 20000, 
             warmup  = 16000,
             chains  = 4,
             cores   = 4,
             inits   = 0,
             control = list(adapt_delta = .99, max_treedepth=13)
)

saveRDS(LV_vero,file="LV_vero.RDS")

LV_trcy<-brm(formula = lotka_volterra,
             prior   = prior,
             data    = trcy_focal,
             family  = poisson("identity"),
             iter    = 20000, 
             warmup  = 16000,
             chains  = 4,
             cores   = 4,
             inits   = 0,
             control = list(adapt_delta = .99, max_treedepth=13)
)



saveRDS(LV_trcy,file="LV_trcy.RDS")


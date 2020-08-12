#beverton holt model fits for competition Vero and Trcy in OPEN environments
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

#vero as focal
 BEV_vero<-brm( formula = beverton_holt,
                prior   = prior,
                data    = vero_focal,
                family  = poisson(link="identity"),
                iter    = 20000,
                warmup  = 16000, 
                cores   = 4,
                chains  = 4,
                inits   = 0 ,
                control =list(adapt_delta = .99, max_treedepth=13)
 )

 saveRDS(BEV_vero,file="BEV_vero.RDS")

#trcy as focal 
BEV_trcy<-brm(formula = beverton_holt,
              prior   = prior,
              data    = trcy_focal,
              family  = poisson(link="identity"),
              iter    = 20000, 
              warmup  = 16000,
              cores   = 4,
              chains  = 4,
              inits   = 0,
              control = list(adapt_delta = .99, max_treedepth=13)
)

saveRDS(BEV_trcy,file="BEV_trcy.RDS")







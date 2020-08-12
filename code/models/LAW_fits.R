library(brms)
source("code/models/set_priors.R")
source("code/clean_data.R")
#the law model
law<- bf(totalseeds~ (exp(lambda)) / ((1 + (alphaii*conspecifics) + (alphaij*heterospecifics))^beta),
          lambda  ~ 1 + env,
          alphaii ~ 1 + env,
          alphaij ~ 1 + env,
          beta    ~ 1 + env,
          nl=T)

#with vero as focal
LAW_vero<-brm(formula = law,
              prior   = prior_exp,
              data    = vero_focal,
              family  = poisson(link="identity"),
              iter    = 20000,
              warmup  = 12000,
              inits   = 0,
              cores   = 4, 
              chains  = 2 ,
              control = list(adapt_delta = .99, max_treedepth=13)
)

saveRDS(LAW_vero,file="LAW_vero.RDS")


#with trcy as focal
      LAW_trcy<-brm(formula = law,
                    prior   = prior_exp,
                    data    = trcy_focal,
                    family  = poisson(link="identity"),
                    iter    = 20000,
                    warmup  = 12000, 
                    cores   = 4,
                    inits   = 0,
                    chains  = 2 ,
                    control = list(adapt_delta = .99, max_treedepth=13)
      )
saveRDS(LAW_trcy,file="LAW_trcy.RDS")







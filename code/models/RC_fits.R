#Ricker models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")


RC_vero<-brm(bf(totalseeds~ (lambda * exp(- (alphaii*conspecifics) - (alphaij*heterospecifics)) ),
                 lambda ~ 1 + env,
                 alphaii ~ 1 + env,
                 alphaij ~ 1 + env,
                 nl=T),
              prior=prior,
              data=vero_focal,family=poisson(link="identity"),
              iter = 4000, warmup = 2000, inits = 1,cores = 4,control =list(adapt_delta = .99)
)

saveRDS(RC_vero,file="RC_vero.RDS")



RC_trcy<-brm(bf(totalseeds~ (lambda *exp( -(alphaii*conspecifics)  - (alphaij*heterospecifics)))  ,
                 lambda ~ 1 + env,
                 alphaii ~ 1 + env,
                 alphaij ~ 1 + env,
                 nl=T),
              prior = prior,
              data=trcy_focal,family=poisson(link="identity"),
              iter = 8000, warmup = 4000, cores = 4 ,chains = 2,control =list(adapt_delta = .99)
)

saveRDS(RC_trcy,file="RC_trcy.RDS")

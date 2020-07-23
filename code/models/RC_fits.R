#Ricker models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")


RC_vero<-brm(bf(totalseeds~ (lambdai * exp(- (alphaii*verodensity) - (alphaij*trcydensity)) ),
                 lambdai ~ 1 + env,
                 alphaii ~ 1 + env,
                 alphaij ~ 1 + env,
                 nl=T),
              prior=prior_vero,
              data=vero_focal,family=poisson(link="identity"),
              iter = 4000, warmup = 2000, cores = 4,inits=1,control =list(adapt_delta = .99)
)

saveRDS(RC_vero,file="RC_vero.RDS")



RC_trcy<-brm(bf(totalseeds~ (lambdaj *exp( -(alphaji*verodensity)  - (alphajj*trcydensity)))  ,
                 lambdaj ~ 1 + env,
                 alphaji ~ 1 + env,
                 alphajj ~ 1 + env,
                 nl=T),
              prior = prior_trcy,
              data=trcy_focal,family=poisson(link="identity"),
              iter = 8000, warmup = 4000, cores = 4, inits=1 ,chains = 2,control =list(adapt_delta = .99)
)

saveRDS(RC_trcy,file="RC_trcy.RDS")

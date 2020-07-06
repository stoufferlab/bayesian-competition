#Ricker models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")


RC_vero<-brm(bf(totalseeds~ (lambdai * exp(- (alphaii*verodensity) - (alphaij*trcydensity)) ) + e*env,
                 lambdai ~ 1 ,
                 alphaii ~ 1 ,
                 alphaij ~ 1 ,
                 e       ~ 1 ,
                 nl=T),
              prior=prior_vero,
              data=vero_focal,family=poisson(),
              iter = 5000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(RC_vero,file="RC_vero.RDS")



RC_trcy<-brm(bf(totalseeds~ (lambdaj *exp( -(alphaji*verodensity)  - (alphajj*trcydensity))) + (e*env) ,
                 lambdaj ~ 1 ,
                 alphaji ~ 1 ,
                 alphajj ~ 1 ,
                 e       ~ 1 ,
                 nl=T),
              prior = prior_trcy,
              data=trcy_focal,family=poisson(),
              iter = 4000, warmup = 2000, cores = 4, chains = 4,inits=0,control =list(adapt_delta = .99)
)

saveRDS(RC_trcy,file="RC_trcy.RDS")

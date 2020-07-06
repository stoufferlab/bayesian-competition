#LV models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")


LV_vero<-brm(bf(totalseeds~ (lambdai *(1 - (alphaii*verodensity) - (alphaij*trcydensity)) ) + e*env,
                 lambdai ~ 1 ,
                 alphaii ~ 1 ,
                 alphaij ~ 1 ,
                 e       ~ 1 ,
                 nl=T),
              prior=prior_vero,
              data=vero_focal,family=poisson(),
              iter = 5000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(LV_vero,file="LV_vero.RDS")



LV_trcy<-brm(bf(totalseeds~ (lambdai *(1 - (alphaji*verodensity) - (alphajj*trcydensity)) ) + e*env,
                lambdai ~ 1 ,
                alphaii ~ 1 ,
                alphaij ~ 1 ,
                e       ~ 1 ,
                nl=T),
             prior=prior_vero,
             data=vero_focal,family=poisson(),
             iter = 5000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(LV_trcy,file="LV_trcy.RDS")

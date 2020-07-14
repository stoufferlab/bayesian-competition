#beverton holt model fits for competition Vero and Trcy in OPEN environments
library(brms)
source("code/models/set_priors.R")



BEV_vero<-brm(bf(totalseeds~ (lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity)) ) ,
                  lambdai ~ 1 + env,
                  alphaii ~ 1 + env,
                  alphaij ~ 1 + env,
                  nl=T),
               prior=prior_vero,
               data=vero_focal,family=poisson(link="identity"),
               iter = 4000, warmup = 2000, cores = 4,chains=4 ,control =list(adapt_delta = .99)
)

saveRDS(BEV_vero,file="BEV_vero.RDS")

BEV_trcy<-brm(bf(totalseeds~ (lambdaj /(1 + (alphaji*verodensity)  + (alphajj*trcydensity)))  ,
                 lambdaj ~ 1 + env,
                 alphaji ~ 1 + env,
                 alphajj ~ 1 + env,
                 nl=T),
              prior = prior_trcy,
              data=trcy_focal,family=poisson(link="identity"),
              iter = 8000, warmup = 4000, cores = 4,chains = 2,control =list(adapt_delta = .99)
)

saveRDS(BEV_trcy,file="BEV_trcy.RDS")







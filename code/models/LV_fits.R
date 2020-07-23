#LV models for VERO vs TRCY
library(brms)
source("code/models/set_priors.R")


# LV_vero<-brm(bf(totalseeds~ lambdai *(1 - alphaii*verodensity - alphaij*trcydensity) ,
#                  lambdai ~ 1 + env ,
#                  alphaii ~ 1 + env,
#                  alphaij ~ 1 + env,
#                  nl=T),
#               prior=prior_vero,
#               data=vero_focal,family=poisson("identity"),
#               iter = 8000, warmup = 4000,chains = 2,cores = 4,control =list(adapt_delta =0.99,max_treedepth=15)
# )
# 
# 
# saveRDS(LV_vero,file="LV_vero.RDS")
# 


LV_trcy<-brm(bf(totalseeds~ lambdaj *(1 - (alphaji*verodensity) - (alphajj*trcydensity)) ,
                lambdaj ~ 1 + env,
                alphajj ~ 1 + env,
                alphaji ~ 1 + env,
                nl=T),
             prior=prior_trcy,
             data=trcy_focal,family=poisson("identity"),
             iter = 8000, warmup = 4000, chains=2,cores = 4,inits = 1,control =list(adapt_delta = .99, max_treedepth=15)
)

saveRDS(LV_trcy,file="LV_trcy.RDS")

library(brms)
source("code/models/set_priors.R")



LAW_vero<-brm(bf(totalseeds~ (lambdai) / ((1 + (alphaii*verodensity) + (alphaij*trcydensity))^bi)  ,
                 lambdai ~ 1 + env,
                 alphaii ~ 1 + env,
                 alphaij ~ 1 + env,
                 bi      ~ 1 + env,
                 nl=T),
              prior=prior_vero_exp,
              data=vero_focal,family=poisson(link="identity"),
              iter = 4000, warmup = 2000, cores = 4,inits = 1, chains=4 ,control =list(adapt_delta = .99)
)

saveRDS(LAW_vero,file="LAW_vero.RDS")

LAW_trcy<-brm(bf(totalseeds~ (lambdaj) /((1 + (alphaji*verodensity)  + (alphajj*trcydensity))^bj ) ,
                 lambdaj ~ 1 + env,
                 alphaji ~ 1 + env,
                 alphajj ~ 1 + env,
                 bj      ~ 1 + env,
                 nl=T),
              prior = prior_trcy,
              data=trcy_focal,family=poisson(link="identity"),
              iter = 4000, warmup = 2000, cores = 4,chains = 4,control =list(adapt_delta = .99)
)

saveRDS(BEV_trcy,file="BEV_trcy.RDS")







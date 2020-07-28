library(brms)
source("code/models/set_priors.R")



LAW_vero<-brm(bf(totalseeds~ lambda / ((1 + (alphaii*conspecifics) + (alphaij*heterospecifics))^b)  ,
                 lambda ~ 1 + env,
                 alphaii ~ 1 + env,
                 alphaij ~ 1 + env,
                 b     ~ 1 + env,
                 nl=T),
              prior=prior_exp,
              data=vero_focal,family=poisson(link="identity"),
              iter = 8000, warmup = 4000, cores = 4, chains=2 ,control =list(adapt_delta = .99)
)

saveRDS(LAW_vero,file="LAW_vero.RDS")


LAW_trcy<-brm(bf(totalseeds~ lambda / ((1 + (alphaii*conspecifics) + (alphaij*heterospecifics))^b)  ,
                 lambda ~ 1 + env,
                 alphaii ~ 1 + env,
                 alphaij ~ 1 + env,
                 b     ~ 1 + env,
                 nl=T),
              prior=prior_exp,
              data=trcy_focal,family=poisson(link="identity"),
              iter = 8000, warmup = 4000, cores = 4, inits = 1,chains=2 ,control =list(adapt_delta = .99)
)
saveRDS(LAW_trcy,file="LAW_trcy.RDS")







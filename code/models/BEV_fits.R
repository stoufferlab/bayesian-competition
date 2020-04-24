#beverton holt model fits for competition Vero and Trcy in OPEN environments
library(brms)
source("code/clean_data.R")
rm(vero_plde)

vero_focal<-vero_trcy[which(vero_trcy$focal=="V"),]
trcy_focal<-vero_trcy[which(vero_trcy$focal=="T"),]

# The multilevel model that I do not want to use because I am not sure it is the right approach

 BEV_vero_multi<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
                       lambdai ~ 1 + (1|place|label),
                       alphaii ~ 1  + (1|place|label),
                       alphaij ~1  + (1|place|label),
                      nl=T),
                    data=vero_focal,family=poisson(),
                    iter = 5000, warmup = 3000, cores = 4, inits=0,control =list(adapt_delta = .99)
 )
 saveRDS(BEV_vero,file="BEV_vero_multi.RDS")



# formula_prior<- bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
#                       lambdai ~ 1 + (1|place|label),
#                       alphaii ~ 1  + (1|place|label),
#                        alphaij ~1  + (1|place|label),
#                       nl=T)
# formula_prior_2<- bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
#                    lambdai ~ 1,
#                    alphaii ~ 1,
#                    alphaij ~1,
#                    nl=T)
# pp<-get_prior(formula = bf(formula_prior_2),family = poisson(),data=vero_focal)
# 
# 
# 

pp_2 = c(
  prior(normal(0,10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 100), nlpar = "lambdai")
)


BEV_vero<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
                  lambdai ~ 1 ,
                  alphaii ~ 1  ,
                  alphaij ~1 ,
                  nl=T),
               prior = pp_2,
               data=vero_focal,family=poisson(),
               iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(BEV_vero,file="BEV_vero.RDS")


pp_3 = c(
  prior(normal(0,10), nlpar = "alphaji"),
  prior(normal(0, 10), nlpar = "alphajj"),
  prior(normal(0, 100), nlpar = "lambdaj")
)

BEV_trcy<-brm(bf(totalseeds~ lambdaj /(1 + (alphaji*verodensity) + (alphajj*trcydensity) ),
                 lambdaj ~ 1 ,
                 alphaji ~ 1  ,
                 alphajj ~1 ,
                 nl=T),
              prior = pp_3,
              data=trcy_focal,family=poisson(),
              iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(BEV_trcy,file="BEV_trcy.RDS")





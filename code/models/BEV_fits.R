#beverton holt model fits for competition Vero and Trcy in OPEN environments
library(brms)
source("code/clean_data.R")
rm(vero_plde)

vero_focal<-vero_trcy[which(vero_trcy$focal=="V"),]
trcy_focal<-vero_trcy[which(vero_trcy$focal=="T"),]

# The multilevel model that I do not want to use because I am not sure it is the right approach

 # BEV_vero_multi<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
 #                       lambdai ~ 1 + (1|place|label),
 #                       alphaii ~ 1  + (1|place|label),
 #                       alphaij ~1  + (1|place|label),
 #                      nl=T),
 #                    data=vero_focal,family=poisson(),
 #                    iter = 5000, warmup = 3000, cores = 4, inits=0,control =list(adapt_delta = .99)
 # )
 # saveRDS(BEV_vero,file="BEV_vero_multi.RDS")



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

prior_vero = c(
  prior(normal(0,10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 100), nlpar = "lambdai")
)


BEV_vero<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
                  lambdai ~ 1 ,
                  alphaii ~ 1  ,
                  alphaij ~1 ,
                  nl=T),
               prior = prior_vero,
               data=vero_focal,family=poisson(),
               iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(BEV_vero,file="BEV_vero.RDS")

prior_trcy = c(
  prior(normal(0,10), nlpar = "alphajj"),
  prior(normal(0, 10), nlpar = "alphaji"),
  prior(normal(0, 10), nlpar = "lambdaj")
)

BEV_trcy<-brm(bf(totalseeds~ lambdaj /(1 + (alphaji*verodensity) + (alphajj*trcydensity) ),
                 lambdaj ~ 1 ,
                 alphaji ~ 1  ,
                 alphajj ~1 ,
                 nl=T),
              prior = prior_trcy,
              data=trcy_focal,family=poisson(),
              iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(BEV_trcy,file="BEV_trcy.RDS")








BEV_trcy_multi<-brm(bf(totalseeds~ lambdaj /(1 + (alphaji*verodensity) + (alphajj*trcydensity) ),
                 lambdaj ~ 1 + (1|place|label) ,
                 alphaji ~ 1 + (1|place|label) ,
                 alphajj ~1  + (1|place|label),
                 nl=T),
              data=trcy_focal,family=poisson(),
              iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(BEV_trcy_multi,file="BEV_trcy_multi.RDS")



library(ggplot2)
ggplot(vero_focal) + geom_line(mapping=aes(x=verodensity,y=totalseeds,col=exposedenv))
ggplot(trcy_focal) + geom_line(mapping=aes(x=trcydensity,y=totalseeds,col=exposedenv))

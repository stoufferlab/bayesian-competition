#beverton holt model fits for competition Vero and Trcy in OPEN environments
library(brms)
source("code/clean_data.R")

vero_focal<-vero_trcy[which(vero_trcy$focal=="V"),]
trcy_focal<-vero_trcy[which(vero_trcy$focal=="T"),]

# 
# BEV_vero<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
#                       lambdai ~ 1 + (1|place|label),
#                       alphaii ~ 1  + (1|place|label),
#                       alphaij ~1  + (1|place|label),
#                       nl=T),
#                    data=vero_focal,family=poisson(),
#                    iter = 5000, warmup = 3000, cores = 4, inits=0,control =list(adapt_delta = .99)
# )
# saveRDS(BEV_vero,file="BEV_vero.RDS")


pp = c(
  prior(normal(0,10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 100), nlpar = "lambdai")
)


BEV_vero_no<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
                  lambdai ~ 1 ,
                  alphaii ~ 1  ,
                  alphaij ~1 ,
                  nl=T),
               prior = pp,
               data=vero_focal,family=poisson(),
               iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)
saveRDS(BEV_vero_no,file="BEV_vero_no.RDS")





vero_open<-vero_trcy[which(vero_trcy$focal=="V" & vero_trcy$exposedenv== "open"),]
vero_woody<-vero_trcy[which(vero_trcy$focal=="V" & vero_trcy$exposedenv== "woody"),]
trcy_open<-vero_trcy[which(vero_trcy$focal=="T" & vero_trcy$exposedenv== "open"),]
trcy_woody<-vero_trcy[which(vero_trcy$focal=="T" & vero_trcy$exposedenv== "woody"),]



BEV_vero_open<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
                lambdai ~ 1 + (1|place|label),
                alphaii ~ 1  + (1|place|label),
                alphaij ~1  + (1|place|label),
                nl=T),
             data=vero_open,family=poisson(),
             iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)
saveRDS(BEV_vero_open,file="BEV_vero_open.RDS")



BEV_vero_woody<-brm(bf(totalseeds~ lambdai /(1 + (alphaii*verodensity) + (alphaij*trcydensity) ),
                      lambdai ~ 1 + (1|place|label),
                      alphaii ~ 1  + (1|place|label),
                      alphaij ~1  + (1|place|label),
                      nl=T),
                   data=vero_woody,family=poisson(),
                   iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)
saveRDS(BEV_vero_woody,file="BEV_vero_woody.RDS")




BEV_trcy_open<-brm(bf(totalseeds~ lambdaj /(1 + (alphaji*verodensity) + (alphajj*trcydensity) ),
                      lambdaj ~ 1 + (1|place|label),
                      alphaji ~ 1  + (1|place|label),
                      alphajj ~1  + (1|place|label),
                      nl=T),
                   data=trcy_open,family=poisson(),
                   iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)
saveRDS(BEV_trcy_open,file="BEV_trcy_open.RDS")



BEV_trcy_woody<-brm(bf(totalseeds~ lambdaj /(1 + (alphaji*verodensity) + (alphajj*trcydensity) ),
                      lambdaj ~ 1 + (1|place|label),
                      alphaji ~ 1  + (1|place|label),
                      alphajj ~1  + (1|place|label),
                      nl=T),
                   data=trcy_woody,family=poisson(),
                   iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)
saveRDS(BEV_trcy_woody,file="BEV_trcy_woody.RDS")




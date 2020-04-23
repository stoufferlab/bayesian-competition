library(brms)
source("code/clean_data.R")

#as vero as the focal, 3 simple models of competition
vero<-vero_trcy[which(vero_trcy$focal=="V"),]
trcy<-vero_trcy[which(vero_trcy$focal=="T"),]
#the beverton hol model
LV_vero<-brm(bf(totalseeds~ lambdai *(1 -(alphaii*verodensity) - (alphaij*trcydensity) ),
                      lambdai ~ 1 + (1|place|label),
                      alphaii ~ 1  + (1|place|label),
                      alphaij ~1  + (1|place|label),
                      nl=T),
                   data=vero,family=poisson(),
                   iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)
saveRDS(LV_vero,file="LV_vero.RDS")



LV_trcy<-brm(bf(totalseeds~ lambdaj *(1 - (alphajj*trcydensity) - (alphaji*verodensity) ),
                      lambdaj ~ 1 + (1|place|label),
                      alphajj ~ 1  + (1|place|label),
                      alphaji ~1  + (1|place|label),
                      nl=T),
                   data=trcy,family=poisson(),
                   iter = 4000, warmup = 2000, cores = 4, inits=0,control =list(adapt_delta = .99)
)

saveRDS(LV_trcy,file="LV_trcy.RDS")
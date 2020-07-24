###get the posterior growth rates for each model fit
require(tidyverse)
require(brms)
source("code/growth_rates.R")
#We have fixed survival and germination for each species. i=vero, j=trcy
gi<-.372
si<-.556
gj<-.258
sj<-.033



#upload the model fits
BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")

RC_vero  <- readRDS("~/bayesian-competition/RC_vero.RDS")
RC_trcy  <- readRDS("~/bayesian-competition/RC_trcy.RDS")


LAW_vero  <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LAW_trcy  <- readRDS("~/bayesian-competition/LAW_trcy.RDS")




#matrix that has a column for each model (the growth rate calculated under that model) and 8000 nrows, which is the size of the posterior distribution for all models
vero_growth<-tibble(.rows = 8000)
trcy_growth<-tibble(.rows = 8000)


#BEVERTON-HOLT growth rates###################################################################
#This could probably be automated but idk 
lambda_vero_BV <- posterior_samples(BEV_vero,pars="b_lambdai_Intercept",exact_match = TRUE)[,1]
growth_vero_BV <- bev_growth(si,gi,lambda_vero_BV)
vero_growth$BV <- growth_vero_BV


rm(lambda_vero_BV,growth_vero_BV)

lambda_trcy_BV  <- posterior_samples(BEV_trcy, pars = "b_lambdaj_Intercept", exact_match = TRUE)[,1]
growth_trcy_BV <- bev_growth(sj,gj,lambda_trcy_BV)
trcy_growth$BV  <- growth_trcy_BV

rm(lambda_trcy_BV,growth_trcy_BV)


###############################################################################
#RICKER model
lambda_vero_RC <- posterior_samples(RC_vero,pars="b_lambdai_Intercept",exact_match = TRUE)[,1]
growth_vero_RC <- ricker_growth(si,gi,lambda_vero_RC)
vero_growth$RC <- growth_vero_RC


rm(lambda_vero_RC,growth_vero_RC)

lambda_trcy_RC  <- posterior_samples(RC_trcy, pars = "b_lambdaj_Intercept", exact_match = TRUE)[,1]
growth_trcy_RC <- ricker_growth(sj,gj,lambda_trcy_RC)
trcy_growth$RC  <- growth_trcy_RC

rm(lambda_trcy_RC,growth_tracy_RC)

#######################################################################################
#LAW model

lambda_vero_LA <- posterior_samples(LAW_vero,pars="b_lambdai_Intercept",exact_match = TRUE)[,1]
b_vero_LA      <- posterior_samples(LAW_vero,pars="b_bi_Intercept",exact_match = TRUE)[,1] 
growth_vero_LA <- law_growth(si,gi,lambda_vero_LA, b_vero_LA)
vero_growth$LA <- growth_vero_LA


rm(lambda_vero_LA,growth_vero_LA)

lambda_trcy_LA <- posterior_samples(LAW_trcy, pars = "b_lambdaj_Intercept", exact_match = TRUE)[,1]
b_trcy_LA      <- posterior_samples(LAW_trcy,pars="b_bj_Intercept",exact_match = TRUE)[,1] 
growth_trcy_LA <- law_growth(sj,gj,lambda_trcy_LA,b_trcy_LA)
trcy_growth$LA  <- growth_trcy_LA

rm(lambda_trcy_LA,growth_trcy_LA)





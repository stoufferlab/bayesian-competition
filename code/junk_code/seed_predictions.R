library(brms)
source("code/models/set_priors.R")

pars<-fixed_model(BEV_vero)

total_seeds<-function(lambda,alphai,alphaj,spi,spj){
  BEV<-lambda/ ( 1 + (alphai*spi) + (alphaj*spj))
  seeds<- exp(BEV)
  return(seeds)
}

vv<-seq(0,30,1)
pred<-total_seeds(pars$lambdai_Intercept,pars$alphaii_Intercept,pars$alphaij_Intercept,vv,0)

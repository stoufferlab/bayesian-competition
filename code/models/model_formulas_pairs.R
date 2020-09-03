###code to fit different models with brms
library(brms)

beverton_holt<-bf(
  totalseeds~ ( exp(lambda) /(1 + (alphaii*conspecifics) + (alphaij*heterospecifics)) ) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)

ricker<- bf(
  totalseeds~ exp(lambda) *exp( -(alphaii*conspecifics)  - (alphaij*heterospecifics))  ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)


hassell<- bf(
  totalseeds~ (exp(lambda)) / ((1 + (alphaii*conspecifics) + (alphaij*heterospecifics) )^  exp(beta)),
         lambda  ~ 1 + env,
         alphaii ~ 1 + env,
         alphaij ~ 1 + env,
         beta    ~ 1 + env,
         nl=T)


lotka_volterra<-bf(
  totalseeds~  exp(lambda) * (1 - (alphaii*conspecifics) - (alphaij*heterospecifics) ),
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)



stouffer <-bf(
  totalseeds~ (exp(lambda) * (1 - ( alphaii*conspecifics) - (alphaij*heterospecifics) ) ) * 
              ( (1- exp( exp(beta)* (conspecifics + 1) )) / ( exp(beta) * (conspecifics+1) )) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  beta    ~ 1 + env,
  nl=T
)




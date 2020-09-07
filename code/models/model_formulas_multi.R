

beverton_holt_multi<-bf(
  totalseeds~ ( exp(lambda) /(1 + (alphaii*conspecifics) + (alphaij*heterospecifics) + (alphaik*totalother)) ) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  alphaik ~ 1 + env,
  nl=T)


ricker_multi<- bf(
  totalseeds~ exp(lambda) *exp( -(alphaii*conspecifics)  - (alphaij*heterospecifics) - (alphaik*totalother))  ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  alphaik ~ 1 + env,
  nl=T)


hassell_multi<- bf(
  totalseeds~ (exp(lambda)) / ((1 + (alphaii*conspecifics) + (alphaij*heterospecifics) + (alphaik*totalother) )^  exp(beta)),
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  alphaik ~ 1 + env,
  beta    ~ 1 + env,
  nl=T)



lotka_volterra_multi<-bf(
  totalseeds~  exp(lambda) * (1 - (alphaii*conspecifics) - (alphaij*heterospecifics) - (alphaik*totalother))  ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaik ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)





stouffer_multi <-bf(
  totalseeds~ (exp(lambda) * (1 - ( alphaii*conspecifics) - (alphaij*heterospecifics) - (alphaik*totalother) ) ) * 
    ( (1- exp( exp(beta)* (conspecifics + 1) )) / ( exp(beta) * (conspecifics+1) )) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  beta    ~ 1 + env,
  alphaik ~ 1 + env,
  nl=T
)

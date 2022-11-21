
# brms nonlinear model structures for various phenomenological models of density-dependent fecundity

require(brms)

# INTERACTION-FREE MODEL
interaction_free<-bf(
  totalseeds~ ( exp(lambda) ) ,
  lambda  ~ 1 + env,
  nl=T)

# BEVERTON-HOLT MODEL
beverton_holt_multi<-bf(
  totalseeds~ ( exp(lambda) /(1 + (alphaii*conspecifics) + (alphaij*heterospecifics) + (alphaik*totalother)) ) ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  alphaik ~ 1 + env,
  nl=T)

# RICKER MODEL
ricker_multi<- bf(
  totalseeds~ exp(lambda) *exp( -(alphaii*conspecifics)  - (alphaij*heterospecifics) - (alphaik*totalother))  ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  alphaik ~ 1 + env,
  nl=T)

# LOTKA-VOLTERRA MODEL
lotka_volterra_multi<-bf(
  totalseeds~  exp(lambda) * (1 - (alphaii*conspecifics) - (alphaij*heterospecifics) - (alphaik*totalother))  ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaik ~ 1 + env,
  alphaij ~ 1 + env,
  nl=T)

# HASSELL MODEL
hassell_multi<- bf(
  totalseeds~ (exp(lambda)) / ((1 + (alphaii*conspecifics) + (alphaij*heterospecifics) + (alphaik*totalother) )^  exp(beta)),
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  alphaik ~ 1 + env,
  beta    ~ 1 + env,
  nl=T)

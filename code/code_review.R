library(brms)
source("code/clean_data.R")

#We will be working with vero as focal so we remove the other species
rm(trcy_focal)


#two models with the same priors
prior_pairs = c(
  prior(normal(0, 10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 10), nlpar = "lambda")
)

#two different functional forms of competition
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

# we fit the models 
beverton_fit <- brm(
  formula = beverton_holt,
  prior = prior_pairs,
  data =  vero_focal,
  family = poisson(link="identity"),
  iter = 4000,
  warmup = 2000,
  cores   = 4,
  chains  = 4,
  inits   = 1,)


ricker_fit <- brm(
  formula = ricker,
  prior = prior_pairs,
  data =  vero_focal,
  family = poisson(link="identity"),
  iter = 4000,
  warmup = 2000,
  cores   = 4,
  chains  = 4,
  inits   = 1,)


### So now I want to compare the models! Which is why I created this ugly function. I want it to be able to give it an undetermined number of models



model_compare<-function(...){
  models <- list(...)
  waic_value <- c()
  looic_value <-c()
  
  for(i in 1:length(models)){
    w <- waic(models[[i]])
    wa <- w$waic
    waic_value[i] <-wa
    
    l <- loo(models[[i]])
    lo <- l$looic
    looic_value[i] <-lo
  }
  
  waic_weights <- model_weights(..., weights = "waic")
  loo_weights <- model_weights(..., weights = "loo")
  out <- cbind(waic_value, waic_weights, looic_value, loo_weights)
  out <- out[order(waic_weights),]
  return(out)
  
}


#We test it, and a bunch of warnings come up, but just as an example

model_compare(beverton_fit,ricker_fit)



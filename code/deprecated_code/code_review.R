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
  nl=TRUE)

ricker<- bf(
  totalseeds~ exp(lambda) *exp( -(alphaii*conspecifics)  - (alphaij*heterospecifics))  ,
  lambda  ~ 1 + env,
  alphaii ~ 1 + env,
  alphaij ~ 1 + env,
  nl=TRUE)

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



model_compare <- function(..., sort.by = "waic"){
  models <- list(...)
  
  # print a warning if only one model is supplied
  if (length(models) == 1) {
    warning("Only one model is supplied. Comparing one model is not meaningful.")
  }
  
  # switching to lapply from for-loop
  ic <- lapply(models, function(m) {
    w <- waic(m)
    wa <- w$estimates["waic", "Estimate"]  # in newer brms access with "estimate", switch back to $waic if you're still using older brms

    l <- loo(m)
    lo <- l$estimates["looic", "Estimate"]

    return(c(waic = wa, looic = lo))
  })
  ic <- do.call(rbind, ic)

  waic_weights <- model_weights(..., weights = "waic")
  loo_weights <- model_weights(..., weights = "loo")
  
  out <-
    cbind(waic = ic[, "waic"],
          waic_weights,
          looic = ic[, "looic"],
          loo_weights)
  out <- out[order(out[, sort.by]),]  # sort.by can be alternatively "looic" now
  return(out)
  
}


#We test it, and a bunch of warnings come up, but just as an example
# both sort methods are the same, since they have the same ranking, but flexible for other scenarios
model_compare(beverton_fit, ricker_fit, sort.by = "waic")
model_compare(beverton_fit, ricker_fit, sort.by = "looic")
model_compare(beverton_fit)

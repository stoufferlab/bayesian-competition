###Model fits for species vero 
library(brms)
source("code/clean_data.R")
#we remove trcy just to be sure
rm(trcy_focal)
source("code/models/set_priors.R")
source("code/models/model_formulas.R")

#Beverton-Holt model
#with pairs
vero_beverton_pairs_zip <- brm(
  formula = beverton_holt_pairs,
  prior = prior_pairs,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)
vero_beverton_pairs_zip<-add_criterion(vero_beverton_pairs_zip, criterion = "waic")
saveRDS(vero_beverton_pairs_zip, file = "model_objects/vero_beverton_pairs_zip.RDS")

#with other species present
vero_beverton_multi_zip <- brm(
  formula = beverton_holt_multi,
  prior = prior_multi,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  save_all_pars = TRUE,
  control = list(adapt_delta = .99, max_treedepth=15)
)

vero_beverton_multi_zip<-add_criterion(vero_beverton_multi_zip, criterion = "waic")
saveRDS(vero_beverton_multi_zip, file = "model_objects/vero_beverton_multi_zip.RDS")


######## The Ricker model
#pairs
vero_ricker_pairs_zip <- brm(
  formula = ricker_pairs,
  prior = prior_pairs,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)
vero_ricker_pairs_zip <- add_criterion(vero_ricker_pairs_zip, criterion = "waic")
saveRDS(vero_ricker_pairs_zip, file ="model_objects/vero_ricker_pairs_zip.RDS")

#with other species present
vero_ricker_multi_zip <- brm(
  formula = ricker_multi,
  prior = prior_multi,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)
vero_ricker_multi_zip <- add_criterion(vero_ricker_multi_zip, criterion = "waic")
saveRDS(vero_ricker_multi_zip, file ="model_objects/vero_ricker_multi_zip.RDS")

####The hassell model

vero_hassel_pairs_zip <- brm(
  formula = hassell_pairs,
  prior = prior_exp_pairs,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)

vero_hassell_pairs_zip <- add_criterion(vero_hassell_pairs_zip, criterion = "waic")
saveRDS(vero_ricker_pairs_zip, file ="model_objects/vero_ricker_pairs_zip.RDS")


###with multi species


vero_hassel_multi_zip <- brm(
  formula = hassell_multi,
  prior = prior_exp_multi,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits = 1,
  control = list(adapt_delta = .99, max_treedepth=15)
)
vero_hassell_multi_zip <- add_criterion(vero_hassell_multi_zip, criterion = "waic")
saveRDS(vero_ricker_multi_zip, file ="model_objects/vero_ricker_multi_zip.RDS")

### lotka volterra


vero_lotka_pairs_zip <- brm(
  formula = lotka_volterra_pairs,
  prior = prior_pairs,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)

vero_lotka_pairs_zip <- add_criterion( vero_lotka_pairs_zip, criterion="waic")
saveRDS(vero_lotka_pairs_zip, file ="model_objects/vero_lotka_pairs_zip.RDS")

#### with multiple species


vero_lotka_multi_zip <- brm(
  formula = lotka_volterra_multi,
  prior = prior_multi,
  data = vero_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)


vero_lotka_multi_zip <- add_criterion( vero_lotka_multi_zip, criterion="waic")
saveRDS(vero_lotka_multi_zip, file ="model_objects/vero_multi_pairs_zip.RDS")






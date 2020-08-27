###Model fits for species trcy 
library(brms)
source("code/clean_data.R")
#we remove trcy just to be sure

source("code/models/set_priors.R")
source("code/models/model_formulas.R")

#Beverton-Holt model
#with pairs
trcy_beverton_pairs_zip <- brm(
  formula = beverton_holt_pairs,
  prior = prior_pairs,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)
trcy_beverton_pairs_zip<-add_criterion(trcy_beverton_pairs_zip, criterion = "waic")
saveRDS(trcy_beverton_pairs_zip, file = "model_objects/trcy_beverton_pairs_zip.RDS")

#with other species present
trcy_beverton_multi_zip <- brm(
  formula = beverton_holt_multi,
  prior = prior_multi,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  save_all_pars = TRUE,
  control = list(adapt_delta = .99, max_treedepth=15)
)

trcy_beverton_multi_zip<-add_criterion(trcy_beverton_multi_zip, criterion = "waic")
saveRDS(trcy_beverton_multi_zip, file = "model_objects/trcy_beverton_multi_zip.RDS")


######## The Ricker model
#pairs
trcy_ricker_pairs_zip <- brm(
  formula = ricker_pairs,
  prior = prior_pairs,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)
trcy_ricker_pairs_zip <- add_criterion(trcy_ricker_pairs_zip, criterion = "waic")
saveRDS(trcy_ricker_pairs_zip, file ="model_objects/trcy_ricker_pairs_zip.RDS")

#with other species present
trcy_ricker_multi_zip <- brm(
  formula = ricker_multi,
  prior = prior_multi,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)
trcy_ricker_multi_zip <- add_criterion(trcy_ricker_multi_zip, criterion = "waic")
saveRDS(trcy_ricker_multi_zip, file ="model_objects/trcy_ricker_multi_zip.RDS")

####The hassell model

trcy_hassel_pairs_zip <- brm(
  formula = hassell_pairs,
  prior = prior_exp_pairs,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)

trcy_hassell_pairs_zip <- add_criterion(trcy_hassell_pairs_zip, criterion = "waic")
saveRDS(trcy_ricker_pairs_zip, file ="model_objects/trcy_ricker_pairs_zip.RDS")


###with multi species


trcy_hassel_multi_zip <- brm(
  formula = hassell_multi,
  prior = prior_exp_multi,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits = 1,
  control = list(adapt_delta = .99, max_treedepth=15)
)
trcy_hassell_multi_zip <- add_criterion(trcy_hassell_multi_zip, criterion = "waic")
saveRDS(trcy_ricker_multi_zip, file ="model_objects/trcy_ricker_multi_zip.RDS")

### lotka volterra


trcy_lotka_pairs_zip <- brm(
  formula = lotka_volterra_pairs,
  prior = prior_pairs,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)

trcy_lotka_pairs_zip <- add_criterion( trcy_lotka_pairs_zip, criterion="waic")
saveRDS(trcy_lotka_pairs_zip, file ="model_objects/trcy_lotka_pairs_zip.RDS")

#### with multiple species


trcy_lotka_multi_zip <- brm(
  formula = lotka_volterra_multi,
  prior = prior_multi,
  data = trcy_focal,
  family = zero_inflated_poisson(link = "identity"),
  iter = 8000,
  warmup = 4000,
  cores   = 4,
  chains  = 4,
  inits   = 1 ,
  control = list(adapt_delta = .99, max_treedepth=15)
)


trcy_lotka_multi_zip <- add_criterion( trcy_lotka_multi_zip, criterion="waic")
saveRDS(trcy_lotka_multi_zip, file ="model_objects/trcy_multi_pairs_zip.RDS")






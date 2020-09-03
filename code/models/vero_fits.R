###Model fits for species vero 
library(brms)
source("code/clean_data.R")
source("code/models/set_priors.R")
source("code/models/model_formulas_pairs.R")
source("code/models/automate_fit.R")

### We create lists of formulas on which to iterate
formulas_pairs <- list( bh = beverton_holt, rc = ricker, hs = hassell, lv = lotka_volterra, st =stouffer)
formulas_multi <- list( bh = beverton_holt_multi, rc = ricker_multi, hs = hassell_multi, lv = lotka_volterra_multi, st =stouffer_multi)

#For species pairs, using a poisson distribution
model_fits(
  focal = "vero",
  data = vero_focal,
  distribution = poisson(link = "identity"),
  formulas = formulas_pairs,
  priors = prior_pairs,
  priors_exponent = prior_exp_pairs,
  num_species = "pairs",
  last_name = "poisson"
)

# For species pairs again, using a zero inflated poisson
model_fits(
  focal = "vero",
  data = vero_focal,
  distribution = zero_inflated_poisson(link = "identity"),
  formulas = formulas_pairs,
  priors = prior_pairs,
  priors_exponent = prior_exp_pairs,
  num_species = "pairs",
  last_name = "zipoisson"
)


#### for multispecies using poisson


model_fits(
  focal = "vero",
  data = vero_focal,
  distribution = poisson(link = "identity"),
  formulas = formulas_multi,
  priors = prior_multi,
  priors_exponent = prior_exp_multi,
  num_species = "multispecies",
  last_name = "poisson"
)



# for multispecies using zip

model_fits(
  focal = "vero",
  data = vero_focal,
  distribution = zero_inflated_poisson(link = "identity"),
  formulas = formulas_multi,
  priors = prior_multi,
  priors_exponent = prior_exp_multi,
  num_species = "multispecies",
  last_name = "poisson"
)




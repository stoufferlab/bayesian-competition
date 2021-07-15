###Model fits for species trcy
library(brms)
source("code/clean_data.R")
source("code/models/set_priors.R")
source("code/models/model_formulas_multi.R")
source("code/models/automate_fit.R")

### We create lists of formulas on which to iterat
formulas <- list(bh = beverton_holt_multi ,rc= ricker_multi)
formulas_lotka <- list(lv = lotka_volterra_multi)

model_fits(
   focal = "trcy",
   data = trcy_focal,
   distribution = poisson(link = "identity"),
   formulas = formulas,
   priors = prior_multi,
   priors_exponent = prior_exp_multi,
   num_species = "multispecies",
   last_name = "poisson"
)

model_fits(
   focal = "trcy",
   data = trcy_focal,
   distribution = poisson(link = "identity"),
   formulas = formulas_lotka,
   priors = prior_lotka,
   priors_exponent = prior_exp_multi,
   num_species = "multispecies",
   last_name = "poisson"
)




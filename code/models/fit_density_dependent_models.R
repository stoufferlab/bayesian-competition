###Model fits for species vero 
require(here)
source(here("code/clean_data.R"))
source(here("code/models/set_priors.R"))
source(here("code/models/model_formulas_multi.R"))
source(here("code/models/automate_fit.R"))

### We create lists of formulas over which to iterate
formulas <- list(
  bh = beverton_holt_multi,
  rc = ricker_multi
)

# fit the models to vero as focal
model_fits(
  focal = "vero",
  data = vero_focal,
  distribution = poisson(link = "identity"),
  formulas = formulas,
  priors = prior_multi,
  num_species = "multispecies",
  last_name = "poisson"
)

# fit the models to trcy as focal
model_fits(
  focal = "trcy",
  data = trcy_focal,
  distribution = poisson(link = "identity"),
  formulas = formulas,
  priors = prior_multi,
  num_species = "multispecies",
  last_name = "poisson"
)

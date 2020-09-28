###Model fits for species vero 
library(brms)
source("code/clean_data.R")
source("code/models/set_priors.R")
source("code/models/model_formulas_multi.R")
source("code/models/automate_fit.R")

### We create lists of formulas on which to iterat
formulas <- list(hs = hassell_multi)


 model_fits(
   focal = "vero",
   data = vero_focal,
   distribution = poisson(link = "identity"),
   formulas = formulas,
   priors = prior_multi,
   priors_exponent = prior_exp_multi,
   num_species = "multispecies",
   last_name = "poisson"
 )




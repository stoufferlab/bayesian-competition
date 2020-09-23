###Model fits for species vero 
library(brms)
source("code/clean_data.R")
source("code/models/set_priors.R")
source("code/models/model_formulas_pairs.R")
source("code/models/automate_fit.R")

### We create lists of formulas on which to iterat
formulas_pairs <- list(bh=beverton_holt,lv=lotka_volterra,rc=ricker,hs = hassell)


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




###Model fits for species trcy 
library(brms)
source("code/clean_data.R")
source("code/models/set_priors.R")

source("code/models/model_formulas_pairs.R")
source("code/models/automate_fit.R")

formulas_pairs <- list(lv=lotka_volterra, hs= hassell)





model_fits(
   focal = "trcy",
   data = trcy_focal,
   distribution = poisson(link = "identity"),
   formulas = formulas_pairs,
   priors = prior_pairs,
   priors_exponent = prior_exp_pairs,
   num_species = "pairs",
   last_name = "poisson"
)







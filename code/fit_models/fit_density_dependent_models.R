
# packages used for parallelization
library(doMC)
registerDoMC(cores = 2)

# for our own functions
require(here)
source(here("code/fit_models/clean_data.R"))
source(here("code/fit_models/brms_model_formulas.R"))
source(here("code/fit_models/brms_model_priors.R"))
source(here("code/fit_models/brms_model_inits.R"))

# we will only fit models listed below
models <- c(
  "beverton_holt",
  "ricker",
  "interaction_free"
)

# iterate over models to be fit
foreach(model_name = models) %:%
  # fit each model for each focal species
  foreach(focal = c("vero", "trcy")) %dopar% {
    # use the appropriate dataset
    if(focal=="vero"){
      data <- vero_focal
    }else{
      data <- trcy_focal
    }

    message(paste("fitting",model_name,"model to focal species",focal,collapse=" "))
    
    # fit the model with brms
    brmsmodel <- brm(
      formula = model_formula(model_name),
      prior = model_prior(model_name),
      family = negbinomial(link = "identity"),
      data =  data,
      # silent=2,
      # refresh=0,
      init = model_init(model_name, data, chains=4),
      iter = 51000,
      warmup = 50000,
      cores   = 4,
      chains  = 4,
      control = list(adapt_delta=0.999, max_treedepth=25)
    )

    # save the model object
    id <- paste0("code/fit_models/model_objects/",focal,"_",model_name,".rds")
    saveRDS(brmsmodel, file = here(id))
  # }
}

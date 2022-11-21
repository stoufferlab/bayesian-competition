
require(here)
source(here("code/clean_data.R"))
source(here("code/models/brms_model_formulas.R"))
source(here("code/models/brms_model_priors.R"))

# we will only fit models listed below
formulas <- list(
  null = interaction_free,
  bevertonholt = beverton_holt_multi,
  ricker = ricker_multi
)

# fit each model for each focal species
for(focal in c("vero", "trcy")){
  if(focal=="vero")
    data <- vero_focal
  else
    data <- trcy_focal

  # iterate over models to be fit
  for(i in seq_along(formulas)){
    # which priors should be used
    if(names(formulas)[i] == 'null')
      prior <- prior_null
    else
      prior <- prior_multi
    
    # fit the model with brms
    model <- brm(
      formula = formulas[[i]],
      prior = prior,
      data =  data,
      family = poisson(link = "identity"),
      iter = 4000,
      warmup = 2000,
      cores   = 4,
      chains  = 4,
      init = 0,
      control = list(adapt_delta = .99, max_treedepth=15)
    )

    # save the model object
    id <- paste0("code/models/model_objects/",focal, "_", names(formulas)[i],".rds")
    saveRDS(model, file = here(id))
  }
}

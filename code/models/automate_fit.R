
#the list of  model formulas to fit 
#function to automate model fits for a single species
# focal is the identifier name of the species, either vero or trcy
# data is the data of the focal, should read from "code/clean_data"
# distribution is the likelihood used to fit the model
# formulast is a list of brms formulas for all the models to fit, with each element having a distinct name
# priors is a list with brms priors for each model parameter,
# priors exp should have an extra parameter for the b parameter in hassell model
# num species is to indicate if you are fitting species pairs or a multispecies system


model_fits <-
  function(focal = "species_name",
           data,
           distribution = poisson(link = "identity"),
           formulas,
           priors,
           priors_exponent,
           num_species = "pairs",
           last_name = "poisson") {
    
    for (i in 1:length(formulas)) {
      if (names(formulas)[[i]] == "hs" | names(formulas)[[i]] == "st") {
        model <- brm(
          formula = formulas[[i]],
          prior = priors_exponent,
          data =  data,
          family = distribution,
          iter = 4000,
          warmup = 2000,
          cores   = 4,
          chains  = 4,
          inits   = 1,
          control = list(adapt_delta = .99, max_treedepth=13)
        )
        
        id <- paste0("model_objects/",focal, "_", names(formulas)[[i]], "_", num_species, "_", last_name,".rds")
        saveRDS(model, file = id)
        
      } else{
        model <- brm(
          formula = formulas[[i]],
          prior = priors,
          data =  data,
          family = distribution,
          iter = 4000,
          warmup = 2000,
          cores   = 4,
          chains  = 4,
          inits   = 1,
          control = list(adapt_delta = .99, max_treedepth=13)
        )
        
        id <- paste0("model_objects/",focal, "_", names(formulas)[[i]], "_", num_species, "_", last_name,".rds")
        saveRDS(model, file = id)
        
      }
      
      
    }
    
  }
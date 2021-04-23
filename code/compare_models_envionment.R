source("code/read_models.R")

library(tidyverse)
library(xtable)
library(brms)

model_compare_environment<- function(..., env){
  models <- list(...)
  
  # print a warning if only one model is supplied
  if (length(models) == 1) {
    warning("Only one model is supplied. Comparing one model is not meaningful.")
  }
  
  #asuming all models were fit with the same data
  data <- models[[1]]$data
  
  #in the control conditions
  
  obs_ind <- which(data$env==env)
 
  
  loos <- loo_subsample(...,
                             observations=obs_ind)
  
   compare <- loos$diffs %>% as.data.frame() %>% select(c(elpd_diff, looic))
  
   compare$weights <- exp(-0.5*compare$elpd_diff)/ sum(exp(-0.5*compare$elpd_diff))

  
  return(compare)
  
}






vero_sunny <- model_compare_environment(vero_bh_multispecies_poisson.rds,
                       vero_lv_multispecies_poisson.rds,
                       vero_rc_multispecies_poisson.rds,
                       env=0) 


vero_woody <- model_compare_environment(vero_bh_multispecies_poisson.rds,
                                        vero_lv_multispecies_poisson.rds,
                                        vero_rc_multispecies_poisson.rds,
                                        env=1)

trcy_sunny <- model_compare_environment(trcy_bh_multispecies_poisson.rds,
                                        trcy_lv_multispecies_poisson.rds,
                                        trcy_rc_multispecies_poisson.rds,
                                        env = 0)

trcy_woody <- model_compare_environment(trcy_bh_multispecies_poisson.rds,
                                        trcy_lv_multispecies_poisson.rds,
                                        trcy_rc_multispecies_poisson.rds,
                                        env = 1)

vero_one <- vero_sunny %>% 
  round(digits = 2) %>%
  select(looic, weights) %>% 
  rownames_to_column(var = "Model")

colnames(vero_one) <-c("Model", "LOOIC sunny", "Weights sunny")



vero_two <- vero_woody %>% 
  round(digits = 2) %>%
  select(looic, weights)%>% 
  rownames_to_column(var = "Model")


colnames(vero_two) <-c("Model", "LOOIC woody", "Weights woody")


vero <-left_join(vero_one, vero_two, by="Model")

vero <- cbind("Species"="Vellia r.", vero)


# trcy --------------------------------------------------------------------


trcy_one <- trcy_sunny %>% 
  round(digits = 2) %>%
  select(looic, weights) %>% 
  rownames_to_column(var = "Model")

colnames(trcy_one) <-c("Model", "LOOIC sunny", "Weights sunny")



trcy_two <- trcy_woody %>% 
  round(digits = 2) %>%
  select(looic, weights)%>% 
  rownames_to_column(var = "Model")


colnames(trcy_two) <-c("Model", "LOOIC woody", "Weights woody")


trcy <-left_join(trcy_one, trcy_two, by="Model")

trcy <- cbind("Species"="Trachymene c.", trcy)


values <- rbind(vero,trcy)


xtable(values)






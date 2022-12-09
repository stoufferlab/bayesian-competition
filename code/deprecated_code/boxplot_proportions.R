source("code/gg_theme.R")
library(tidyverse)
library(xtable)

#funciton that for every model combination in our results, calculates the proportion of feasible points in 10 samples, doing it n_samples times, to have boxplots instead of bargraphs when talking about proportions in coexistence
extract_proportion_samples <- function(mod){
  
  mod <- mod %>% unite("combo",vero_model:trcy_model )
  
  model_combos <- unique(mod$combo) %>% as.list()
  
 proportions <- lapply(model_combos, function(x, 
                                              n_samples){
   one_combo <- mod %>% filter(combo == x)
   
   # nn <- nrow(one_combo)
   # feasible_samples <- c()
   # for( i in 1:n_samples){
   #   sample <- one_combo[sample(1:nn, 10),]
   #   feasible_proportion <- mean(sample$feasibility)
   # #  print(feasible_proportion)
   #   feasible_samples <- c(feasible_samples, feasible_proportion)
   # } 
   
   prop<- mean( one_combo$feasibility) %>% round(digits = 2)

   results <- data.frame("model"= x,
                         "proportion"= prop)
   return(results)
   
 },n_samples=n_samples) 
 
  all_models <- do.call(rbind, proportions)
  
  all_models  <- all_models %>% separate(model, into=c("vero_model", "trcy_model"),sep = "_")
  
  return(all_models)
}


results_sunny_bounded <- readRDS("~/bayesian-competition/integration_objects/results_sunny_bounded.RDS")
results_woody_bounded <- readRDS("~/bayesian-competition/integration_objects/results_woody_bounded.RDS")


sunny <- extract_proportion_samples(mod = results_sunny_bounded)
woody <- extract_proportion_samples(mod = results_woody_bounded)

tab <- cbind(sunny, woody$proportion) 



colnames(tab)<- c("Vellia r. model", "Trachymene c. model", "Coexistence proportion", "Coexistence proportion woody")


xtable(tab)

source("code/gg_theme.R")
library(tidyverse)


#funciton that for every model combination in our results, calculates the proportion of feasible points in 10 samples, doing it n_samples times, to have boxplots instead of bargraphs when talking about proportions in coexistence
extract_proportion_samples <- function(mod,
                                       n_samples){
  
  mod <- mod %>% unite("combo",vero_model:trcy_model )
  
  model_combos <- unique(mod$combo) %>% as.list()
  
 proportions <- lapply(model_combos, function(x, 
                                              n_samples){
   one_combo <- mod %>% filter(combo == x)
   nn <- nrow(one_combo)
   feasible_samples <- c()
   for( i in 1:n_samples){
     sample <- one_combo[sample(1:nn, 10),]
     feasible_proportion <- mean(sample$feasibility)
   #  print(feasible_proportion)
     feasible_samples <- c(feasible_samples, feasible_proportion)
   } 
   
   results <- data.frame("model"= x,
                         "proportions"= feasible_samples,
                         "mean_proportion"= mean(feasible_samples))
   return(results)
   
 },n_samples=n_samples) 
 
  all_models <- do.call(rbind, proportions)
  
  all_models  <- all_models %>% separate(model, into=c("vero_model", "trcy_model"),sep = "_")
  
  return(all_models)
}



pp <- extract_proportion_samples(mod = results_sunny_bounded,n_samples = 100)
pp <- pp %>% mutate(value=1)


pp<-pp[which(pp$vero_model== "Lotka-Volterra"| pp$vero_model=="Ricker"),]

ggplot(pp) +
  geom_boxplot(mapping = aes(x=value,
                             y=proportions),
               fill="mediumseagreen") +
  theme_alba +
  facet_grid(vero_model~trcy_model)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(breaks = NULL, labels = NULL)+
  labs(x="")

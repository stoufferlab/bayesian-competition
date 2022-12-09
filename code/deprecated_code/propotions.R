library(ggplot2)
library(ggpubr)

# figures for the paper ---------------------------------------------------

source("code/access_results.R")

sunny_results <- filter(final_results, environment=="open") %>% filter( vero_model != "Lotka-Volterra")%>% filter( trcy_model != "Lotka-Volterra" )
woody_results <- filter(final_results, environment=="woody")  %>% filter( vero_model != "Lotka-Volterra")%>% filter( trcy_model != "Lotka-Volterra" )

proportions_model <- function(mod){
  mod <- mod %>% unite("combo", vero_model:trcy_model, remove = FALSE)
  combos <- unique(mod$combo)
  
  props <- lapply(seq_along(combos), function(x){
    one_combo <- mod %>% filter(combo == combos[x])
  
     props_one <-sum(one_combo$feasibility)/nrow(one_combo) 
     
     results <- data.frame("vero_model"= unique(one_combo$vero_model),
                           "trcy_model"= unique(one_combo$trcy_model),
                           "proportion"= round(props_one,digits = 2)
                             )
     return(results)
    
  })
  
  final <- do.call(rbind, props)
  return(final)
}


proportions_model(sunny_results)

proportions_model(woody_results)

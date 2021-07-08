source("code/compare_models_envionment.R")
source("code/boxplot_proportions.R")



values$Model <- c("Beverton-Holt",
                  "Ricker",
                  "Lotka-Volterra",
                  "Beverton-Holt",
                  "Lotka-Volterra",
                  "Ricker")

colnames(values)<-c("species","model","looic1","sunny","looic2" ,"woody")

colnames(tab)<- c("Vellia r.", "Trachymene c.", "sunny", "woody")


final<- c()

for(i in 1:nrow(tab)){
  one_row <- tab[i,]
  
  sunny_weight_vellia <- values %>% filter(species== "Vellia r.") %>%
    filter( model == one_row$`Vellia r.`) %>% select( sunny)
  
  sunny_weight_trcy <- values %>% filter(species=="Trachymene c.") %>%
    filter( model == one_row$`Trachymene c.`) %>% select( sunny)
  
  
  woody_weight_vellia <- values %>% filter(species== "Vellia r.") %>%
    filter( model == one_row$`Vellia r.`) %>% select(woody)
  
  woody_weight_trcy <- values %>% filter(species=="Trachymene c.") %>%
    filter( model == one_row$`Trachymene c.`) %>% select( woody)
  
  sunny_prop <- one_row$sunny %>% as.numeric() 
  woody_prop <- one_row$woody %>% as.numeric()
  

  
  s_p <- sunny_weight_vellia * sunny_weight_trcy * sunny_prop  
  w_p <- woody_weight_vellia * woody_weight_trcy * woody_prop
  
  
  colnames(s_p) <- NULL
  colnames(w_p) <- NULL
  
  results <- data.frame("Vellia r." = one_row$`Vellia r.`,
                        "Trachymene c."= one_row$`Trachymene c.`,
                        "proportion_sunny"= sunny_prop,
                        "probability_sunny"= s_p,
                        "proportion_woody"= woody_prop,
                        "probability_woody"= w_p)
  
  final <- rbind(final, results)
  
}



table3 <- final[,3:ncol(final)] %>% round(digits = 2)


tt <-cbind( final[1:2], table3)


xtable::xtable(tt)

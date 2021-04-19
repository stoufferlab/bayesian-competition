source("code/read_models.R")
source("code/compare_models.R")

library(tidyverse)
library(xtable)


vero <- model_compare(vero_bh_multispecies_poisson.rds, 
                      vero_lv_multispecies_poisson.rds, 
                      vero_rc_multispecies_poisson.rds) %>%
  round(digits = 3) %>% 
  as.data.frame() 

vero <- cbind(Species = "Vellia rosea", Model=c("Beverton-Holt", "Ricker", "Lotka-Volterra"),vero)

rownames(vero) <- NULL
colnames(vero)<- c("Species", "Model","WAIC", "WAIC weights", "LOOIC", "LOOIC weights")


# trcy --------------------------------------------------------------------

trcy <- model_compare(trcy_bh_multispecies_poisson.rds, 
                      trcy_lv_multispecies_poisson.rds, 
                      trcy_rc_multispecies_poisson.rds) %>%
  round(digits = 4) %>% 
  as.data.frame() 

trcy <- cbind(Species = "Trachymene cyanopetala", Model=c("Beverton-Holt",  "Lotka-Volterra", "Ricker"),trcy)


rownames(trcy) <- NULL
colnames(trcy)<- c("Species", "Model","WAIC", "WAIC weights", "LOOIC", "LOOIC weights")
models <- rbind(vero, trcy)
rownames(models) <- NULL


xtable(models)


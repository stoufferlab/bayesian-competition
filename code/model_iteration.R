#iterating over models

library(brms)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#We source everything known to human kind...


source("code/read_models.R")
source("code/feasibility_toolbox.R")
source("code/saavedra_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")
source("code/model_combo.R")

#survival and germination for Vero (i) and Trcy(j)
#gi<-.372
#si<-.556
#gj<-.258
#sj<-.033


#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324


# the list of models over to iterate
vero_models <- list( vero_bh_multispecies_poisson.rds,
                    # vero_lv_multispecies_poisson.rds,
                     vero_rc_multispecies_poisson.rds)

trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                  # trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds)

model_grid <- combined_models(vero_models = vero_models,
                                   trcy_models = trcy_models,
                                   si =si,
                                   gi =gi,
                                   gj =gj,
                                   sj=sj,
                                   Ni_max = 1e3,
                                   Nj_max =1e3,
                                   env=FALSE,
                                   bounded = TRUE)

saveRDS(model_grid,
        file = "results_sunny_bounded.RDS")


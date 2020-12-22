#iterating over models

library(brms)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#We source everything known to human kind...

source("code/gg_theme.R")
source("code/read_models.R")
source("code/model_toolbox.R")
source("code/model_combo.R")
source("code/integration_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")


#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033

# the list of models over to iterate
vero_models <- list( vero_bh_multispecies_poisson.rds,
                     vero_lv_multispecies_poisson.rds,
                     vero_rc_multispecies_poisson.rds)

trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                   trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds)

model_grid_sunny<- combined_models(vero_models = vero_models,
                                   trcy_models = trcy_models,
                                   si =si,
                                   gi =gi,
                                   gj =gj,
                                   sj=sj,
                                   Ni = 1e4,
                                   Nj =1e4,
                                   env=FALSE,
                                   make_plot = TRUE)

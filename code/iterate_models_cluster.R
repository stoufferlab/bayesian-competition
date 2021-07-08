
# cluster script ----------------------------------------------------------

library(brms)
library(tidyverse)
library(rlist)

#We source everything known to human kind...

source("code/read_models.R")
source("code/feasibility_toolbox.R")
source("code/saavedra_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")


# Upload models and relevant params --------------------------------------------------------

#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324

Ni_max = 5e3
Nj_max =5e3
env=FALSE
bounded = TRUE

# the list of models over to iterate
vero_models <- list( vero_bh_multispecies_poisson.rds,
                     vero_lv_multispecies_poisson.rds,
                     vero_rc_multispecies_poisson.rds)

trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                   trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds)


# Define iterat ---------------------------------------------------------


args <- commandArgs(trailingOnly=T)
i <- as.numeric(args[1])


vero_models_names <- lapply(vero_models, function(x){return(x$name)}) %>% unlist()
trcy_models_names <- lapply(trcy_models, function(x){return(x$name)}) %>% unlist()
  
model_names_combo <- expand.grid(vero_models_names,
                                   trcy_models_names)

vero <- list.filter(vero_models, name==model_names_combo[i,1] )[[1]]
trcy <- list.filter(trcy_models, name==model_names_combo[i,2] )[[1]]
    
  posterior_feasibility(vero_model=vero,
                          trcy_model = trcy,
                          si = si,
                          gi = gi,
                          sj = sj,
                          gj = gj,
                          Ni_max  = Ni_max,
                          Nj_max = Nj_max,
                          env = env,
                          bounded = bounded)
    
  
  
  


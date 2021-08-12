

require("brms")
require("tidyverse")

#We source everything known to human kind...

source("code/read_models.R")
source("code/feasibility_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")
source("code/parameter_sensitivity_functions.R")


# Upload models and relevant params --------------------------------------------------------

#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324

Ni_max = 5e3
Nj_max =5e3


# BH and BH ---------------------------------------------------------------

bh_bh_open <-
  posterior_feasibility_sensitivity(
    vero_model = vero_bh_multispecies_poisson.rds,
    trcy_model = trcy_bh_multispecies_poisson.rds,
    si = si,
    gi = gi,
    Ni_max = Ni_max,
    sj = sj,
    gj = gj,
    Nj_max = Nj_max,
    env = FALSE
  )


bh_bh_woody <-
  posterior_feasibility_sensitivity(
    vero_model = vero_bh_multispecies_poisson.rds,
    trcy_model = trcy_bh_multispecies_poisson.rds,
    si = si,
    gi = gi,
    Ni_max = Ni_max,
    sj = sj,
    gj = gj,
    Nj_max = Nj_max,
    env = TRUE
  )


bh_bh_sensitivity <- rbind(bh_bh_open,
                           bh_bh_woody)

write_csv(bh_bh_sensitivity, file = "bh_bh_sensitiviyt.csv")


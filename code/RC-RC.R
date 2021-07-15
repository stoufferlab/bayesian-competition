
# my own way of parallel ----------------------------------------------------------

require("brms")
require("tidyverse")
require("rlist")
require("sp")
require("SpatialGraph")
require("tidyverse")
require("polylabelr")
require("tidyverse")
require("mvtnorm")
require("cubature")


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
bounded = TRUE



# Define iterat ---------------------------------------------------------




posterior_feasibility(vero_model=vero_rc_multispecies_poisson.rds,
                      trcy_model = trcy_rc_multispecies_poisson.rds,
                      si = si,
                      gi = gi,
                      sj = sj,
                      gj = gj,
                      Ni_max  = Ni_max,
                      Nj_max = Nj_max,
                      env = FALSE,
                      bounded = bounded)





posterior_feasibility(vero_model=vero_rc_multispecies_poisson.rds,
                      trcy_model = trcy_rc_multispecies_poisson.rds,
                      si = si,
                      gi = gi,
                      sj = sj,
                      gj = gj,
                      Ni_max  = Ni_max,
                      Nj_max = Nj_max,
                      env = TRUE ,
                      bounded = bounded)




source("code/read_models.R")
source("code/feasibility_toolbox.R")
source("code/saavedra_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")

# 
# gi<-.372
# si<-.556
# gj<-.258
# sj<-.033


#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324

ptm <- proc.time()

 test_2<- posterior_feasibility(vero_model = vero_bh_multispecies_poisson.rds,
                              trcy_model = trcy_bh_multispecies_poisson.rds,
                              si = si,
                              gi = gi,
                              sj = sj,
                              gj = gj,
                              Ni_max  = 1e3,
                              Nj_max = 1e3,
                              env = FALSE,
                              bounded = TRUE)



proc.time() - ptm
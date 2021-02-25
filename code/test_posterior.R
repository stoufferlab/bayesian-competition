source("code/read_models.R")
source("code/integration_toolbox.R")
source("code/model_toolbox.R")
source("code/determine_radius.R")
source("code/determine_boundary.R")

gi<-.372
si<-.556
gj<-.258
sj<-.033

ptm <- proc.time()

 test_2<- posterior_feasibility(vero_model = vero_lv_multispecies_poisson.rds,
                              trcy_model = trcy_bh_multispecies_poisson.rds,
                              si = si,
                              gi = gi,
                              sj = sj,
                              gj = gj,
                              Ni = 1e4,
                              Nj = 1e4,
                              env = FALSE,
                              make_plot = TRUE)



proc.time() - ptm
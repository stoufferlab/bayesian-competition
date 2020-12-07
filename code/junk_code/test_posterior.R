#quick test to see how the posterior distributions behave under the new code. Only using 500 draws per model and arbritrarly large Ns

#testing how the code works with mean parameter values
source("code/contrained_toolbox.R")
source("code/model_toolbox.R")
source("code/read_models.R")
source("code/model_combo.R")

#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033

vero_models <- list( vero_bh_multispecies_poisson.rds,
                     vero_lv_multispecies_poisson.rds,
                     vero_rc_multispecies_poisson.rds)

trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                   trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds)



post<- combined_models( vero_models = vero_models,
  trcy_models = trcy_models,
                      si = si,
                      gi = gi,
                      Ni = 1000,
                      sj = sj,
                      gj = gj,
                      Nj = 1000,
                      env = FALSE,
                      make_plot = FALSE)

saveRDS(post, file="posterior_feasibility_test.RDS")



#testing how the code works with mean parameter values
source("code/contrained_toolbox.R")
source("code/model_toolbox.R")
source("code/read_models.R")

#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033




trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                   trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds)

#Using abundances of 1000 for both species, this is how using the mean paramter values work
#not automated because its just for kicks


par(mfrow=c(3,3))



lapply(trcy_models,function(m){
  posterior_feasibility( vero_model=vero_bh_multispecies_poisson.rds,
                         trcy_model = m,
                         si = si,
                         gi = gi,
                         Ni = 1000,
                         sj = sj,
                         gj = gj,
                         Nj = 1000,
                         env = FALSE )
}
      )




lapply(trcy_models,function(m){
  posterior_feasibility( vero_model=vero_lv_multispecies_poisson.rds,
                         trcy_model = m,
                         si = si,
                         gi = gi,
                         Ni = 1000,
                         sj = sj,
                         gj = gj,
                         Nj = 1000,
                         env = FALSE )
}
)





lapply(trcy_models,function(m){
  posterior_feasibility( vero_model=vero_rc_multispecies_poisson.rds,
                         trcy_model = m,
                         si = si,
                         gi = gi,
                         Ni = 1000,
                         sj = sj,
                         gj = gj,
                         Nj = 1000,
                         env = FALSE )
}
)

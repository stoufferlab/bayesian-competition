library(brms)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#We source everything known to human kind...

source("code/gg_theme.R")
source("code/read_models.R")
source("code/model_toolbox.R")
source("code/model_combo.R")
source("code/feasibility_toolbox.R")


#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033

# the list of models over to iterate
vero_models <- list( vero_bh_multispecies_poisson.rds,
                     vero_lv_multispecies_poisson.rds,
                     vero_rc_multispecies_poisson.rds,
                     vero_hs_multispecies_poisson.rds)

trcy_models<- list(trcy_bh_multispecies_poisson.rds,
                   trcy_lv_multispecies_poisson.rds,
                   trcy_rc_multispecies_poisson.rds,
                   trcy_hs_multispecies_poisson.rds)

model_grid_sunny<- combined_models(vero_models = vero_models,
                                   trcy_models = trcy_models,
                                   si =si,
                                   gi =gi,
                                   gj =gj,
                                   sj=sj,
                                   env=FALSE)



closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }


sample_grid<-function(model_grid){
  omegas<- seq(0,1,0.001)
  thetas<- seq(0,90, .1)
  feasibility_grid<-matrix(0, nrow=length(thetas), ncol = length(omegas))
  colnames(feasibility_grid) <-omegas
  rownames(feasibility_grid) <-thetas
  
  for(i in 1:nrow(model_grid)){
    omega_i <- closest(omegas, model_grid$omega_results[i]) %>% as.character()
    theta_i <- closest(thetas, model_grid$theta_results[i]) %>% as.character()
    #print(omega_i)
    feasibility_grid[theta_i,omega_i] <- feasibility_grid[theta_i,omega_i] +1
  }
  
 return(feasibility_grid)
}


one_grid<- sample_grid(model_grid_sunny)


one_grid_vectors<-one_grid  %>% as.data.frame() %>%
  rownames_to_column("theta") %>%
  pivot_longer(-c(theta), names_to = "omega")



quantile_range <- quantile(one_grid_vectors$value, probs = seq(0, 1, 0.01))
                   
one_grid_vectors$quantile<-findInterval(one_grid_vectors$value, quantile_range, all.inside = TRUE)

one_grid_vectors$theta <- one_grid_vectors$theta %>% as.numeric()
one_grid_vectors$omega <- one_grid_vectors$omega %>% as.numeric()
ggplot(one_grid_vectors)+
  geom_tile(aes(x=omega, y=theta, fill=value))+
  theme_alba+
  scale_fill_continuous(low="grey", high="black")+
  geom_abline(slope = 180, intercept = 0, col="grey50")+
  xlim(0,0.3)


source("code/clean_data.R")
source("code/model_assesment.R")
require(brms)

BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LV_vero <- readRDS("~/bayesian-competition/LV_vero.RDS")

BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")
RC_trcy <- readRDS("~/bayesian-competition/RC_trcy.RDS")
LAW_trcy <- readRDS("~/bayesian-competition/LAW_trcy.RDS")
LV_trcy <- readRDS("~/bayesian-competition/LV_trcy.RDS")

vero_names <-c("results/beverton_vero.pdf", "results/ricker_vero.pdf", "results/law_vero.pdf","results/lotka_vero.pdf")
trcy_names <-c("results/beverton_trcy.pdf", "results/ricker_trcy.pdf", "results/law_trcy.pdf","results/lotka_trcy.pdf")

view_models<-function(m1,m2,m3,m4,data, vector_names){
  models<-list(m1,m2,m3,m4)
  for (i in 1:length(models)){
    pdf( file = vector_names[i],  width = 7, onefile = TRUE)
    m <- models[[i]]
    print(class(m))
   convergence(model = m)
   prediction_quantiles(m,data = data)
   dev.off()
  }
  
}

view_models(BEV_vero, RC_vero, LAW_vero, LV_vero, data= vero_focal, vector_names = vero_names)

view_models(BEV_trcy, RC_trcy, LAW_trcy, LV_trcy, data= trcy_focal, vector_names = trcy_names)



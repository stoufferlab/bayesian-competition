# figures for presentation

source("code/theta_figure.R")

BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")

RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
RC_trcy <- readRDS("~/bayesian-competition/RC_trcy.RDS")

LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LAW_trcy <- readRDS("~/bayesian-competition/LAW_trcy.RDS")

LV_vero <- readRDS("~/bayesian-competition/LV_vero.RDS")
LV_trcy <- readRDS("~/bayesian-competition/LV_trcy.RDS")

#BEVERTON HOLT - BEVERTON HOLT 

#
BEV_BEV_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0, title = "Beverton-Holt & Beverton-Holt")


##  RICKER - BEVERTON HOLT

RC_BEV_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0,  title = "Ricker & Beverton-Holt")


#BEVERTON  HOLT - RICKER
BEV_RC_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0,  title = "Beverton-Holt & Ricker")

##  RICKER - RICKER
RC_RC_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0,  title = "Ricker & Ricker")

##### LAW - BEVERTON HOLT
LAW_BEV_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0,  title = "Law & Beverton-Holt")



#### LAW - RICKER

LAW_RC_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0,  title = "Law & Ricker")

##### BEVERTON -LAw

BEV_LAW_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0,  title = "Beverton-Holt  & Law")


##  RICKER - LAW

RC_LAW_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0,  title = "Ricker & Law")


########### LAW - LAW
LAW_LAW_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0,  title = "Law & Law")

############## LOTKA - BEV

LV_BEV_0 <-make_figure_theta(vero_model =  LV_vero, vero_function = lotka_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0,  title = "Lotka & Beverton-Holt")

#### LOTKA-RC

LV_RC_0 <-make_figure_theta(vero_model =  LV_vero, vero_function = lotka_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0,  title = "Lotka & Ricker")



#### LOTKA -LAW 
LV_LAW_0 <-make_figure_theta(vero_model =  LV_vero, vero_function = lotka_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0,  title = "Lotka & Law")


### LOTKA -LOTKA
LV_LV_0 <-make_figure_theta(vero_model =  LV_vero, vero_function = lotka_growth, vero_exp = 0 , trcy_model = LV_trcy, trcy_function = lotka_growth, trcy_exp = 0, env = 0,  title = "Lotka & Lotka")

###BEVERTON - LOTKA
BEV_LV_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = LV_trcy, trcy_function = lotka_growth, trcy_exp = 0, env = 0,  title = "Beverton-Holt & Lotka")

#### RICKER - LOTKA
RC_LV_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = LV_trcy, trcy_function = lotka_growth, trcy_exp = 0, env = 0,  title = "Ricker & Lotka")

##### LAW - LOTKA

LAW_LV_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = LV_trcy, trcy_function = lotka_growth, trcy_exp = 0, env = 0,  title = "Law & Lotka")



 pdf(file="omega.pdf", width = 14)
 grid.arrange(BEV_BEV_0, RC_BEV_0, LAW_BEV_0, LV_BEV_0,
              BEV_RC_0,  RC_RC_0,  LAW_RC_0,  LV_RC_0,
              BEV_LAW_0, RC_LAW_0, LAW_LAW_0, LV_LAW_0,
              BEV_LV_0,  RC_LV_0,  LAW_LV_0, LV_LV_0,
              ncol=4, nrow=4)
 dev.off()
 
 
 
 

 
 
 
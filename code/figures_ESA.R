# figures for presentation

source("code/theta_figure.R")

BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")

RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
RC_trcy <- readRDS("~/bayesian-competition/RC_trcy.RDS")

LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LAW_trcy <- readRDS("~/bayesian-competition/LAW_trcy.RDS")

#BEVERTON HOLT - BEVERTON HOLT 

#
BEV_BEV_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0)



BEV_BEV_1 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 1)



##  RICKER - BEVERTON HOLT

RC_BEV_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0)

RC_BEV_1 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 1)



#BEVERTON  HOLT - RICKER
BEV_RC_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0)

BEV_RC_1 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 1)


##  RICKER - RICKER

RC_RC_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0)

RC_RC_1 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 1)

##### LAW - BEVERTON HOLT
LAW_BEV_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0)

LAW_BEV_1 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 1)



#### LAW - RICKER

LAW_RC_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 0)

LAW_RC_1 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = RC_trcy, trcy_function = ricker_growth, trcy_exp = 0, env = 1)



##### BEVERTON -LAw

BEV_LAW_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0)


BEV_LAW_1 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 1)


##  RICKER - LAW

RC_LAW_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0)



RC_LAW_1 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 1)


########### LAW - LAW
LAW_LAW_0 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 0)

LAW_LAW_1 <-make_figure_theta(vero_model =  LAW_vero, vero_function = law_growth, vero_exp = 1 , trcy_model = LAW_trcy, trcy_function = law_growth, trcy_exp = 1, env = 1)



 pdf(file="omega_multi_model.pdf", width = 14)
 grid.arrange(BEV_BEV_0, RC_BEV_0, LAW_BEV_0,
              BEV_RC_0,  RC_RC_0,  LAW_RC_0,
              BEV_LAW_0, RC_LAW_0, LAW_LAW_0, 
              ncol=3, nrow=3)
 dev.off()
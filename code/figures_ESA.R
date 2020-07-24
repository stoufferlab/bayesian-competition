# figures for presentation

source("code/theta_figure.R")

BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")

RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")

#BEVERTON HOLT - BEVERTON HOLT 

#
BEV_BEV_0 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0)


BEV_BEV_1 <-make_figure_theta(vero_model =  BEV_vero, vero_function = bev_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 1)



##  RICKER - BEVERTON HOLT

RC_BEV_0 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 0)

RC_BEV_1 <-make_figure_theta(vero_model =  RC_vero, vero_function = ricker_growth, vero_exp = 0 , trcy_model = BEV_trcy, trcy_function = bev_growth, trcy_exp = 0, env = 1)










# 
# pdf(file="omega_more", width = 14)
# grid.arrange(BV_BV, RC_BV, LA_BV,
#              BV_RC, RC_RC, LA_RC,
#              BV_LA, RC_LA, LA_LA, 
#              ncol=3, nrow=3)
# dev.off()
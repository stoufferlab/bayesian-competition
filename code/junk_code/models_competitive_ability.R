library(brms)
library(ggplot2)
library(tidyverse)
library(ggsci)
source("code/competitive_ability.R")

posterior_competitive_ability<-function(model, fun, s ,g, exp_param){
  post        <- posterior_samples(model)
  lambda      <- exp(post$b_lambda_Intercept)
  lambda_env  <- exp(lambda + post$b_lambda_env)
  alphaii     <- post$b_alphaii_Intercept
  alphaii_env <- post$b_alphaii_Intercept + post$b_alphaii_env
  alphaij     <- post$b_alphaij_Intercept
  alphaij_env <- post$b_alphaij_Intercept + post$b_alphaij_env
 
  
  if(exp_param){
    b_param     <- exp(post$b_beta_Intercept)
    b_param_env <- exp(post$b_beta_Intercept + post$b_beta_env)
    ability     <- fun(s,g,lambda,b_param,alphaii,alphaij)
    ability_env <- fun(s,g,lambda_env, b_param, alphaii_env, alphaij_env)
  }else{
    ability     <- fun(s,g,lambda,alphaii,alphaij)
    ability_env <- fun(s,g,lambda_env, alphaii_env, alphaij_env)
  }
  
  all_posterior <- cbind(post, ability, ability_env)
  return(all_posterior)
}



BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LV_vero <- readRDS("~/bayesian-competition/LV_vero.RDS")

gi<-.372
si<-.556
gj<-.258
sj<-.033

bv_ca_vero <- posterior_competitive_ability(model = BEV_vero, 
                                            fun =  bev_competitive_ability, 
                                            s=si,
                                            g=gi,
                                            exp_param=0 ) %>% mutate ( model = "beverton_holt", species= "vero")
lv_ca_vero <- posterior_competitive_ability(model =LV_vero, 
                                            fun =lotka_competitive_ability, 
                                            s=si,
                                            g=gi,
                                            exp_param= 0 ) %>% mutate ( model = "lotka_volterra", species = "vero")
rc_ca_vero <- posterior_competitive_ability(model =RC_vero, 
                                            fun =ricker_competitive_ability, 
                                            s=si,
                                            g=gi,
                                            exp_param= 0 ) %>% mutate ( model = "ricker", species = "vero")
lw_ca_vero <- posterior_competitive_ability(model =LAW_vero, 
                                            fun =law_competitive_ability, 
                                            s=si,
                                            g=gi,
                                            exp_param= 1 )   %>% mutate ( model = "law", species = "vero")


vero_competitive_posterior <- bind_rows(bv_ca_vero, lv_ca_vero, rc_ca_vero, lw_ca_vero)




BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")
RC_trcy <- readRDS("~/bayesian-competition/RC_trcy.RDS")
LAW_trcy <- readRDS("~/bayesian-competition/LAW_trcy.RDS")
LV_trcy <- readRDS("~/bayesian-competition/LV_trcy.RDS")


bv_ca_trcy <- posterior_competitive_ability(model = BEV_trcy, 
                                            fun =  bev_competitive_ability, 
                                            s=sj,
                                            g=gj,
                                            exp_param=0 ) %>% mutate ( model = "beverton_holt", species ="trcy")
lv_ca_trcy <- posterior_competitive_ability(model =LV_trcy, 
                                            fun =lotka_competitive_ability, 
                                            s=sj,
                                            g=gj,
                                            exp_param= 0 ) %>% mutate ( model = "lotka_volterra", species = "trcy")
rc_ca_trcy <- posterior_competitive_ability(model =RC_trcy, 
                                            fun =ricker_competitive_ability, 
                                            s=sj,
                                            g=gj,
                                            exp_param= 0 ) %>% mutate ( model = "ricker", species = "trcy")
lw_ca_trcy <- posterior_competitive_ability(model =LAW_trcy, 
                                            fun =law_competitive_ability, 
                                            s=sj,
                                            g=gj,
                                            exp_param= 1 )   %>% mutate ( model = "law", species = "trcy")



trcy_competitive_posterior <- bind_rows(bv_ca_trcy, lv_ca_trcy, rc_ca_trcy, lw_ca_trcy)


all_competitive_posterior <- bind_rows(vero_competitive_posterior, trcy_competitive_posterior)



ggplot(all_competitive_posterior) + geom_density(aes(x=ability, fill = model), alpha= 0.7 )+ 
  facet_wrap( ~ species) +
  theme_bw() +
  scale_fill_locuszoom() +
  xlim(0,.2) +
  ylim(0,200) +
  labs (x="Sensitivity to competition")
 
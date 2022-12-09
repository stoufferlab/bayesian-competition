###########examining the growth rates 
library(ggplot2)
library(tidyverse)


source("code/model_toolbox.R")
source("code/growth_rates.R")

BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")

RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
RC_trcy <- readRDS("~/bayesian-competition/RC_trcy.RDS")

LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LAW_trcy <- readRDS("~/bayesian-competition/LAW_trcy.RDS")


gi<-.372
si<-.556
gj<-.258
sj<-.033




vero_bev  <- posterior_parameters(model = BEV_vero, fun = bev_growth, s=si, g=gi, exp_param = FALSE) %>%select(growth) %>% mutate ( model = "bev") 

vero_rc   <- posterior_parameters(model = RC_vero, fun = ricker_growth, s=si, g=gi, exp_param = FALSE) %>%select(growth) %>% mutate ( model = "rc") 
vero_law  <- posterior_parameters(model = LAW_vero, fun = law_growth, s=si, g=gi, exp_param = TRUE)  %>%select(growth) %>% mutate ( model = "law") 
vero_law<- vero_law %>% filter (growth < 25000)


growth_vero<-rbind(vero_bev,vero_rc)

c1<-rethinking::col.alpha("seagreen4", .7)
c2<-rethinking::col.alpha("seagreen1", .7)
c3<-rethinking::col.alpha("seagreen3", .7) 

p1<-ggplot(growth_vero, aes(x=growth, fill=model)) + geom_density() + theme_bw() + scale_fill_manual(values=c(c1,c2 ))
 
p2<-ggplot(vero_law, aes(x=growth)) + geom_density(fill=c3) +theme_bw() 

















  
library(brms)
library(ggplot2)
library(ggpubr)

#We source everything known to human kind...

source("code/gg_theme.R")
source("code/read_models.R")
source("code/model_toolbox.R")


#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033



#A function that returns a data frame with posterior equilibrium values ofr each models

multiple_equilibriums <-function(..., s,g) {
   
    models <- list(...)
    posteriors <- lapply(models, function(m) {
      one_posterior <- posterior_parameters(m,s,g)
      return(one_posterior)
    })
    all_posteriors <- do.call(rbind, posteriors)
    
  }

#order is important in this function, each model has its own growth function and equilibrium function, and name


#This is the order : beverton holt (bh), lotka volterra (lv), ricker (rc) and hassell (hs)
model_names <- c("Beverton-Holt", "Lotka-Volterra", "Ricker", "Hassell")


vero <- multiple_equilibriums(
  vero_bh_multispecies_poisson.rds,
  vero_lv_multispecies_poisson.rds,
  vero_rc_multispecies_poisson.rds,
  vero_hs_multispecies_poisson.rds,
  s = si,
  g = gi
)


trcy <- multiple_equilibriums(
  trcy_bh_multispecies_poisson.rds,
  trcy_lv_multispecies_poisson.rds,
  trcy_rc_multispecies_poisson.rds,
  trcy_hs_multispecies_poisson.rds,
  s = sj,
  g = gj
)

trcy$model<-factor(trcy$model, levels = model_names)
vero$model<-factor(vero$model, levels= model_names)


  
vero_0<-ggplot(vero) +
  geom_density(mapping = aes(x = equilibrium, fill = model, linetype=model) ,
               alpha = 0.8, show.legend = FALSE) +
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_fill_manual(values=palette_alba)+
  theme_alba+
  xlim(0,400) +
  ylim(0,0.027)



trcy_0<-ggplot(trcy) +
  geom_density(mapping = aes(x = equilibrium, fill = model, linetype=model) ,
               alpha = 0.8, show.legend = FALSE) +
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_fill_manual(values=palette_alba)+
  theme_alba +
  xlim(0,400)+
  ylim(0,0.027)



vero_1<-ggplot(vero) +
  geom_density(mapping = aes(x = env_equilibrium, fill = model, linetype=model) ,
               alpha = 0.8,show.legend = FALSE) +
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_fill_manual(values=palette_alba) +
  theme_alba +
  xlim(0,1500)+
  ylim(0,0.0125)


trcy_1<-ggplot(trcy) +
  geom_density(mapping = aes(x = env_equilibrium, fill = model, linetype=model) ,
               alpha = 0.8, show.legend = TRUE) +
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_fill_manual(values=palette_alba)+
  theme_alba + 
  xlim(0,1500) +
  ylim(0,0.0125)







all_figures<-ggarrange( vero_0, trcy_0 ,
                        vero_1, trcy_1,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

# 
# env_figures<-ggarrange( vero_1, trcy_1 ,
#                         labels = c("A", "B"),
#                         ncol = 1, nrow = 2)
# 


setwd("/home/alba/bayesian_competiton_ms")
pdf(file = "equilibrium.pdf", width = 7, height = 7/1.6)
annotate_figure(all_figures,
             
                bottom = text_grob("Monoculture equilibrium abundance",
                                  size = 14),
                left = text_grob("Density",  rot = 90, size = 14),
                
)
 
dev.off()

setwd("/home/alba/bayesian-competition")

# 
# 
# pdf(file = "results/env_equilibrium.pdf", width = 6, height = 7.5/1.6)
# annotate_figure(env_figures,
#                 
#                 bottom = text_grob("Monoculture equilibrium abundance",
#                                    size = 14),
#                 left = text_grob("Density",  rot = 90, size = 14),
#                 
# )
# 
# dev.off()
library(brms)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
source("code/gg_theme.R")
source("code/per_capita_loss.R")
source("code/read_models.R")
source("code/model_loss.R")

model_names <- c("Beverton-Holt", "Ricker")
loss_functions <- list(bh_loss, rc_loss)

vero_models <-list(vero_bh_multispecies_poisson.rds,
                   vero_rc_multispecies_poisson.rds)

trcy_models<-list(trcy_bh_multispecies_poisson.rds,
                  trcy_rc_multispecies_poisson.rds)

vero_loss <- multiple_loss(models = vero_models,
                           model_names = model_names,
                           loss_fucntions = loss_functions,
                           neighbors = 60  )

trcy_loss <- multiple_loss(models = trcy_models,
                           model_names = model_names,
                           loss_fucntions = loss_functions,
                           neighbors =  60 )

#we need to tell it the order of stuff, always as the model names
trcy_loss$model<-factor(trcy_loss$model, levels = model_names)
vero_loss$model<-factor(vero_loss$model, levels= model_names)

#trcy_loss <- trcy_loss %>% melt( id.vars = c("N", "model"), value.name = "response")

yy <- 0.1
vero_cons_0 <- ggplot(data=vero_loss, aes(x=N, y=conspecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy)+
  ggtitle("Conspecifics")


vero_het_0 <- ggplot(data=vero_loss, aes(x=N, y=heterospecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = TRUE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy) +
  ggtitle("Heterospecifics")

vero_cons_1 <- ggplot(data=vero_loss, aes(x=N, y=conspecific_response_env)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy)


vero_het_1 <- ggplot(data=vero_loss, aes(x=N, y=heterospecific_response_env)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy)



vero_figures<-ggarrange( vero_cons_0 , vero_het_0, vero_cons_1, vero_het_1,
                        labels = c("A", "B", "C","D"),
                        ncol = 2, nrow = 2)





pdf(file = "../bayesian_competition_ms//vero_loss.pdf", width = 7, height = 7/1.6)
annotate_figure(vero_figures,
                
                bottom = text_grob("Neighbor density ",
                                   size = 12),
                left = text_grob("Per capita seed loss",  rot = 90, size = 12))


dev.off()


trcy_cons_0 <- ggplot(data=trcy_loss, aes(x=N, y=conspecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy)+
  ggtitle("Conspecifics")


trcy_het_0 <- ggplot(data=trcy_loss, aes(x=N, y=heterospecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = TRUE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy) +
  ggtitle("Heterospecifics")

trcy_cons_1 <- ggplot(data=trcy_loss, aes(x=N, y=conspecific_response_env)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy)


trcy_het_1 <- ggplot(data=trcy_loss, aes(x=N, y=heterospecific_response_env)) + 
  geom_line(aes(linetype = model, color= model), size=1, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,yy)



trcy_figures<-ggarrange( trcy_cons_0 , trcy_het_0, trcy_cons_1, trcy_het_1,
                         labels = c("A", "B", "C","D"),
                         ncol = 2, nrow = 2)




pdf(file = "../bayesian_competition_ms//trcy_loss.pdf", width = 7, height = 7/1.6)
annotate_figure(trcy_figures,
                
                bottom = text_grob("Neighbor density ",
                                   size = 12),
                left = text_grob("Per capita loss",  rot = 90, size = 12))


dev.off()

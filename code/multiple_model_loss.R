library(brms)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
source("code/gg_theme.R")
source("code/per_capita_loss.R")
source("code/read_models.R")
source("code/model_loss.R")

model_names <- c("Beverton-Holt", "Lotka-Volterra", "Ricker", "Hassell")
loss_functions <- list(bh_loss, lv_loss, rc_loss, hs_loss)

vero_models <-list(vero_bh_pairs_poisson.rds,
                   vero_lv_pairs_poisson.rds,
                   vero_rc_pairs_poisson.rds,
                   vero_hs_pairs_poisson.rds)

trcy_models<-list(trcy_bh_pairs_poisson.rds,
                  trcy_lv_pairs_poisson.rds,
                  trcy_rc_pairs_poisson.rds,
                  trcy_hs_pairs_poisson.rds)

vero_loss <- multiple_loss(models = vero_models,
                           model_names = model_names,
                           loss_fucntions = loss_functions,
                           neighbors = 200  )

trcy_loss <- multiple_loss(models = trcy_models,
                           model_names = model_names,
                           loss_fucntions = loss_functions,
                           neighbors = 200  )

#we need to tell it the order of stuff, always as the model names
trcy_loss$model<-factor(trcy_loss$model, levels = model_names)
vero_loss$model<-factor(vero_loss$model, levels= model_names)

#trcy_loss <- trcy_loss %>% melt( id.vars = c("N", "model"), value.name = "response")

vero_0 <- ggplot(data=vero_loss, aes(x=N, y=conspecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=.7, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,0.1)


vero_1 <- ggplot(data=vero_loss, aes(x=N, y=heterospecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=.7, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,0.1)





trcy_1 <- ggplot(data=trcy_loss, aes(x=N, y=heterospecific_response)) + 
          geom_line(aes(linetype = model, color= model), size=.7, show.legend = TRUE)+
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,0.1)

trcy_0 <- ggplot(data=trcy_loss, aes(x=N, y=conspecific_response)) + 
  geom_line(aes(linetype = model, color= model), size=.7, show.legend = FALSE)+
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_color_manual(values=palette_alba) +
  theme_alba +
  ylim(0,0.1)
  






all_figures<-ggarrange( vero_0, vero_1, trcy_0, trcy_1 ,
                        labels = c("A", "B", "C","D"),
                        ncol = 2, nrow = 2)




pdf(file = "results/loss.pdf", width = 8, height = 8.5/1.6)
annotate_figure(all_figures,
                
                bottom = text_grob("Neighbor density ",
                                   size = 14),
                left = text_grob("Per capita loss",  rot = 90, size = 14),
                top = text_grob("Conspecifics                             Heterospecifics", size = 12)
                
                
)

dev.off()



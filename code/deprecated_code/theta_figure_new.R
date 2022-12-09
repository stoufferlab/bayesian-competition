library(brms)
library(ggplot2)
library(ggpubr)

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


model_grid_woody<- combined_models(vero_models = vero_models,
                                   trcy_models = trcy_models,
                                   si =si,
                                   gi =gi,
                                   gj =gj,
                                   sj=sj,
                                   env=TRUE)




model_names <- c("Beverton-Holt", "Lotka-Volterra", "Ricker", "Hassell")

model_grid_sunny$vero_model <- factor(model_grid_sunny$vero_model, levels = model_names)
model_grid_sunny$trcy_model <- factor(model_grid_sunny$trcy_model, levels = model_names)


model_grid_woody$vero_model <- factor(model_grid_woody$vero_model, levels = model_names)
model_grid_woody$trcy_model <- factor(model_grid_woody$trcy_model, levels = model_names)



col1 <- rethinking::col.alpha("mediumseagreen", .05)
col2 <- rethinking::col.alpha("grey50", .05)

feasibility_plot_sunny<-ggplot(model_grid_sunny) +
  geom_point(
    mapping = aes(
      x = omega_results,
      y = theta_results,
      col = as.factor(feasibility_results)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = Theta_mean
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  geom_abline(slope = 180, intercept = 0, col="grey50")+
  xlim(0,0.3)+
  ylim(0, 60)+
  facet_grid(trcy_model~vero_model)



feasibility_plot_woody<-ggplot(model_grid_woody) +
  geom_point(
    mapping = aes(
      x = omega_results,
      y = theta_results,
      col = as.factor(feasibility_results)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = Theta_mean
  ), col= "goldenrod1")+
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  xlim(0,0.3)+
  ylim(0,60)+
  facet_grid(trcy_model~vero_model)+
  geom_abline(slope = 180, intercept = 0, col="grey50")


#we save it where the manuscript lives
setwd("/home/alba/bayesian_competiton_ms")


pdf(file="feasibility_sunny_models.pdf", width = 8, height = 8, onefile = FALSE)
annotate_figure(feasibility_plot_sunny,
                
                bottom = text_grob(expression(Omega),
                                   size = 14),
                left = text_grob(expression(theta),  rot = 90, size = 14),
                
)

dev.off()


pdf(file="feasibility_woody_models.pdf", width = 8, height = 8, onefile = FALSE)
annotate_figure(feasibility_plot_woody,
                
                bottom = text_grob(expression(Omega),
                                   size = 14),
                left = text_grob(expression(theta),  rot = 90, size = 14),
                
)

dev.off()


setwd("/home/alba/bayesian-competition")
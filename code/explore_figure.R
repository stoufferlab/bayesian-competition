library(ggplot2)
library(tidyverse)
source("code/gg_theme.R")

mod<-  readRDS("~/bayesian-competition/results_sunny_22jan21.RDS")
unbounded_mod <- readRDS("~/bayesian-competition/results_sunny_unbounded_22jan21.RDS")

col1 <- rethinking::col.alpha("mediumseagreen", .7)
col2 <- rethinking::col.alpha("grey50",0.7)


integration <- ggplot(mod) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = distance,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = distance_mean
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_grid(trcy_model~vero_model)+
   xlim(0,0.13)+
   ylim(-1.5,3.5)


fig1<- annotate_figure(integration,
                       
                       bottom = text_grob(expression(Omega),
                                          size = 12),
                       left = text_grob("Distance",  rot = 90, size = 12),
                       top = text_grob("Velleia rosea",
                                          size = 12,
                                      face = "italic" ),
                       right = text_grob("Trachymene cyanopetala",
                                       size = 12,
                                       face = "italic",
                                       rot = 270),
                       
)





unbounded_integration <- ggplot(unbounded_mod) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = distance,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = distance_mean
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_grid(trcy_model~vero_model)



fig2<- annotate_figure(unbounded_integration,
                       
                       bottom = text_grob(expression(Omega),
                                          size = 12),
                       left = text_grob("Distance",  rot = 90, size = 12),
                       top = text_grob("Velleia rosea",
                                       size = 12,
                                       face = "italic" ),
                       right = text_grob("Trachymene cyanopetala",
                                         size = 12,
                                         face = "italic",
                                         rot = 270),
                       
)


get_proportions <- function(mod) {
  models <- unique(mod$vero_model)
  results <- c()
  
  for(i in 1:length(models)){
    mod1 <- mod[which(mod$vero_model == models[i]),]
    for(j in 1:length(models)){
      mod2 <-mod1[which(mod1$trcy_model==models[j]),]
      print(sum(mod2$feasibility))
      coexistence <- sum(mod2$feasibility)/nrow(mod2)
      name <- paste0(unique(mod2$vero_model)
                     ,"-", 
                     unique(mod2$trcy_model))
      
      row_results <- data.frame("propotion_coexistence" = coexistence,
                                "models" = name)
      results <- rbind(results,
                       row_results)
      
    }
  }
  return(results)
}







proportions_coexistence <- get_proportions(mod =mod)
proportions_coexistence$models_ab <- c("BH-BH", "BH-LV",
                                       "BH-RC", 
                                       "LV-BH",
                                       "LV-LV",
                                       "LV-RC",
                                       "RC-BH",
                                       "RC-LV",
                                       "RC-RC")


proportions_coexistence_unbounded <- get_proportions(mod = unbounded_mod)
proportions_coexistence_unbounded$models_ab <- c("BH-BH", "BH-LV",
                                                 "BH-RC", 
                                                 "LV-BH",
                                                 "LV-LV",
                                                 "LV-RC",
                                                 "RC-BH",
                                                 "RC-LV",
                                                 "RC-RC")

bars <- ggplot(proportions_coexistence)+
  geom_bar(aes(x=models_ab, 
               y= propotion_coexistence),
           stat = "identity",
           fill = col1)+
  theme(axis.text.x = element_text(angle = 0),
        legend.text = element_text(size = 12),
        legend.position = c(0.850, 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = NA),
        panel.border = element_rect(fill = NA,
                                    colour = "grey20"),
        panel.grid = element_line(colour = "white"),
        panel.grid.minor = element_line(size = rel(0.5)),
        strip.background = element_rect(fill = "white",
                                        colour = "grey20"),
        legend.key = element_rect(fill = "white",
                                  colour = NA),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black")) +
  geom_abline(intercept = max(proportions_coexistence$propotion_coexistence),
              slope = 0,
              col="grey50",
              size =1,
              linetype = "dashed")





bars_unbounded <- ggplot(proportions_coexistence_unbounded)+
  geom_bar(aes(x=models_ab, 
               y= propotion_coexistence),
           stat = "identity",
           fill = col1)+
  theme(axis.text.x = element_text(angle = 0),
        legend.text = element_text(size = 12),
        legend.position = c(0.850, 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = NA),
        panel.border = element_rect(fill = NA,
                                    colour = "grey20"),
        panel.grid = element_line(colour = "white"),
        panel.grid.minor = element_line(size = rel(0.5)),
        strip.background = element_rect(fill = "white",
                                        colour = "grey20"),
        legend.key = element_rect(fill = "white",
                                  colour = NA),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black")) +
  geom_abline(intercept = max(proportions_coexistence$propotion_coexistence),
              slope = 0,
              col="grey50",
              size =1,
              linetype = "dashed")



fig1b<- annotate_figure(bars,
                        
                        
                        left = text_grob("Proportion of coexistence",  rot = 90, size = 12),)




fig2b<- annotate_figure(bars_unbounded,
                        
                        
                        left = text_grob("Proportion of coexistence",  rot = 90, size = 12),)



pdf("woody_results.pdf", width = 10, height = 15)
a<- ggarrange(fig1,
              fig1b,
              ncol=1,
              nrow=2,
              heights = c(1,0.5))
a
dev.off()




pdf("sunny_results_unbounded.pdf", width = 10, height = 15)
b<- ggarrange(fig2,
              fig2b,
              ncol=1,
              nrow=2,
              heights = c(1,0.5))
b
dev.off()





ff <-   annotate_figure(a,
                        
                        
                        top = text_grob("Bounded Integration",  size = 12, hjust = 2, face = "bold"),)



ff1 <-   annotate_figure(b,
                         
                         
                         top = text_grob("Unbounded Integration",  size = 12, hjust = 2, face = "bold"),)


pdf("compare_integration.pdf", width = 15, height = 10)
ggarrange(ff,
          ff1,
          ncol=2,
          nrow=1)
# heights = c(1,0.5))

dev.off()







# comparing omegas --------------------------------------------------------



get_combos <- function(mod) {
  models <- unique(mod$vero_model)
  results <- c()
  
  for(i in 1:length(models)){
    mod1 <- mod[which(mod$vero_model == models[i]),]
    for(j in 1:length(models)){
      mod2 <-mod1[which(mod1$trcy_model==models[j]),]
      name <- paste0(unique(mod2$vero_model)
                     ,"-", 
                     unique(mod2$trcy_model))
      
      mod2$name <- name
      results <- rbind(results,
                       mod2)
      
    }
  }
  return(results)
}



mod1 <- get_combos(mod)
mod2 <- get_combos(unbounded_mod)




uno <- ggplot(mod1)+
  geom_point(aes(x = Omega_saaveda,
                 y = Omega,
                 col= name)) +
  geom_abline(intercept = 0, 
              slope = 1,
              col="grey50")+
  scale_color_brewer(palette = "BrBG") +
  xlim(0,0.3)+
  ylim(0,0.3)+
  theme(#legend.title = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 10),
    panel.background = element_rect(fill = "white",
                                    colour = NA),
    panel.border = element_rect(fill = NA,
                                colour = "grey20"),
    panel.grid = element_line(colour = "white"),
    panel.grid.minor = element_line(size = rel(0.5)),
    strip.background = element_rect(fill = "white",
                                    colour = "grey20"),
    legend.key = element_rect(fill = "white",
                              colour = NA),
    strip.text.x = element_text(size = 12, colour = "black"),
    strip.text.y = element_text(size = 12, colour = "black"))



dos<-ggplot(mod2)+
  geom_point(aes(x = Omega_saaveda,
                 y = Omega,
                 col= name)) +
  scale_color_brewer(palette = "BrBG")+
  geom_abline(intercept = 0, 
              slope = 1,
              col="grey50")+ 
  theme(#legend.title = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 10),
    panel.background = element_rect(fill = "white",
                                    colour = NA),
    panel.border = element_rect(fill = NA,
                                colour = "grey20"),
    panel.grid = element_line(colour = "white"),
    panel.grid.minor = element_line(size = rel(0.5)),
    strip.background = element_rect(fill = "white",
                                    colour = "grey20"),
    legend.key = element_rect(fill = "white",
                              colour = NA),
    strip.text.x = element_text(size = 12, colour = "black"),
    strip.text.y = element_text(size = 12, colour = "black"))





tres <- ggplot(mod1)+
  geom_point(aes(x = theta_saavedra,
                 y = distance,
                 col= name)) +
  scale_color_brewer(palette = "BrBG")+
  theme(#legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 10),
    panel.background = element_rect(fill = "white",
                                    colour = NA),
    panel.border = element_rect(fill = NA,
                                colour = "grey20"),
    panel.grid = element_line(colour = "white"),
    panel.grid.minor = element_line(size = rel(0.5)),
    strip.background = element_rect(fill = "white",
                                    colour = "grey20"),
    legend.key = element_rect(fill = "white",
                              colour = NA),
    strip.text.x = element_text(size = 12, colour = "black"),
    strip.text.y = element_text(size = 12, colour = "black"))+
  ylim(-2,4)





cuatro <- ggplot(mod2)+
  geom_point(aes(x = theta_saavedra,
                 y = distance,
                 col= name)) +
  scale_color_brewer(palette = "BrBG")+
  theme(#legend.title = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 10),
    panel.background = element_rect(fill = "white",
                                    colour = NA),
    panel.border = element_rect(fill = NA,
                                colour = "grey20"),
    panel.grid = element_line(colour = "white"),
    panel.grid.minor = element_line(size = rel(0.5)),
    strip.background = element_rect(fill = "white",
                                    colour = "grey20"),
    legend.key = element_rect(fill = "white",
                              colour = NA),
    strip.text.x = element_text(size = 12, colour = "black"),
    strip.text.y = element_text(size = 12, colour = "black"))+
  ylim(-2,4)




uno_a<- annotate_figure(uno,
                        
                        bottom = text_grob(expression(Omega),
                                           size = 12),
                        left = text_grob("Bounded Integration",  rot = 90, size = 12) )



dos_a<- annotate_figure(dos,
                        
                        bottom = text_grob(expression(Omega),
                                           size = 12),
                        left = text_grob("Unbounded Integration",  rot = 90, size = 12) )






tres_a<- annotate_figure(tres,
                         
                         bottom = text_grob(expression(theta),
                                            size = 12),
                         left = text_grob("Distance from edge, bounded",  rot = 90, size = 12) )



cuatro_a<- annotate_figure(cuatro,
                           
                           bottom = text_grob(expression(theta),
                                              size = 12),
                           left = text_grob("Distance from edge, unbounded",  rot = 90, size = 12) )





pdf("compare_params.pdf", width = 15, height = 10)
ggarrange(uno_a,
          dos_a,
          tres_a,
          cuatro_a,
          ncol = 2,
          nrow = 2)
dev.off()





p1 <- ggplot(mod2)+
  geom_point(aes(x = theta_saavedra,
                 y = distance,
                 col= Omega))+
  scale_color_gradient(low = "firebrick1", high = "darkblue")


p2 <- ggplot(mod1)+
  geom_point(aes(x = theta_saavedra,
                 y = distance,
                 col= Omega))+
  scale_color_gradient(low = "firebrick1", high = "darkblue")+
  ylim(-2,4)

col3 <- rethinking::col.alpha("firebrick1", .3)
col4 <- rethinking::col.alpha("darkblue",0.3)

ggplot(mod2)+
  geom_point(aes(x = Omega_saaveda,
                 y = theta_saavedra,
                 col= distance))+
  scale_color_gradient(low = "firebrick1", high = "darkblue")


ggplot(mod1)+
  geom_point(aes(x = Omega_saaveda,
                 y = theta_saavedra,
                 col= distance))+
  scale_color_gradient(low = "firebrick1", high = "darkblue")



ggarrange(p1,
          p2,
          ncol = 2,
          nrow = 1)




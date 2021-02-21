library(ggplot2)
library(tidyverse)
source("code/gg_theme.R")

mod<-  readRDS("~/bayesian-competition/results_sunny_22jan21.RDS")
unbounded_mod <- readRDS("~/bayesian-competition/results_sunny_unbounded_22jan21.RDS")

mod <- mod %>% mutate (distance_norm = (distance_center -distance_growth)) %>% 
  mutate( distance_mean_norm = ( distance_mean_center - distance_mean_growth ))

col1 <- rethinking::col.alpha("mediumseagreen", .3)
col2 <- rethinking::col.alpha("grey50",.9)


integration <- ggplot(mod) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = distance_growth,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = distance_mean_growth
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_grid(trcy_model~vero_model) 

fig1<- annotate_figure(integration,
                       
                       bottom = text_grob(expression(Omega),
                                          size = 12),
                       left = text_grob("Distance (absolute)",  rot = 90, size = 12),
                       top = text_grob("Velleia rosea",
                                          size = 12,
                                      face = "italic" )
                       # ,
                       # right = text_grob("Trachymene cyanopetala",
                       #                 size = 12,
                       #                 face = "italic",
                       #                 rot = 270),
                       
)


integration <- ggplot(mod) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = distance_norm,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = distance_mean_norm
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_grid(trcy_model~vero_model) 

fig2<- annotate_figure(integration,
                         
                         bottom = text_grob(expression(Omega),
                                            size = 12),
                         left = text_grob("Distance (relative)",  rot = 90, size = 12),
                         top = text_grob("Velleia rosea",
                                         size = 12,
                                         face = "italic" ),
                         right = text_grob("Trachymene cyanopetala",
                                           size = 12,
                                           face = "italic",
                                           rot = 270),
                         
)




ggarrange( fig1, fig2, ncol = 2,nrow = 1)




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





dos <- ggplot(mod1)+
  geom_point(aes(x = Omega,
                 y = distance_norm,
                 col= as.factor(feasibility))) +
  #scale_color_brewer(palette = "BrBG") +
  scale_color_manual(values = c(col2, col1))+
  theme(#legend.title = element_text(size = 12),
    legend.text = element_text(size = 8),
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


uno <- annotate_figure(uno,
                
                bottom = text_grob(expression(Omega),
                                   size = 12),
                left = text_grob("Normalized distance",  rot = 90, size = 12)
                
)


dos <- annotate_figure(dos,
                       
                       bottom = text_grob(expression(Omega),
                                          size = 12),
                       left = text_grob("Normalized distance",  rot = 90, size = 12)
                       
)



ggarrange( uno, dos, ncol = 2, nrow = 1)
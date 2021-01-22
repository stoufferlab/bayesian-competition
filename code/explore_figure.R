library(ggplot2)
library(tidyverse)
source("code/gg_theme.R")

mod<- readRDS("~/bayesian-competition/Beverton-HoltandBeverton-Holt.RDS")


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
  xlim(0,0.13)+
  ylim(-1.5,3.5)+
  facet_grid(trcy_model~vero_model)


fig1<- annotate_figure(integration,
                       
                       bottom = text_grob(expression(Omega),
                                          size = 12),
                       left = text_grob("Distance",  rot = 90, size = 12),
                       top = text_grob("Velleia Rosea",
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



fig2<- annotate_figure(bars,
                       
                       
                       left = text_grob("Proportion of coexistence",  rot = 90, size = 12),)


pdf("preeliminary_results.pdf", width = 10, height = 15)
ggarrange(fig1,
          fig2,
          ncol=1,
          nrow=2,
          heights = c(1,0.5))
dev.off()
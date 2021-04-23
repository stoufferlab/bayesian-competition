library(ggplot2)
library(tidyverse)
source("code/gg_theme.R")



make_figure<-function(mod,
                      name){
  

  col1 <- rethinking::col.alpha("mediumseagreen", .3)
  col2 <- rethinking::col.alpha("grey60",.3)
  
  col3<- rethinking::col.alpha("#e3004e", 1)

  
  
  
  integration <- ggplot(mod) +
    geom_point(
      mapping = aes(
        x = proportion,
        y = distance,
        col = as.factor(feasibility)
      ),
      show.legend = FALSE
    ) +
    geom_point( mapping = aes(
      x = proportion_mean,
      y = distance_mean,
      shape= as.factor(feasibility_mean)),
      size=3,
      col=col3,
      show.legend = FALSE ) +
    theme_alba +
    scale_color_manual(values = c(col2, col1)) +
    facet_grid(trcy_model~vero_model)  +
    geom_abline(intercept = 0, slope = 0, linetype ="dashed", col="grey50")+
    xlim(0,1.1)
  
 
  
  fig1<- annotate_figure(integration,
                         
                         bottom = text_grob( "Relative coexistence area",
                                             size = 12),
                      
                         left = text_grob("Distance from the edge, D",
                                          rot = 90, size = 12),
                         top = text_grob("Velleia rosea",
                                         size = 12,
                                         face = "italic" )
                         ,
                         right = text_grob("Trachymene cyanopetala",
                                           size = 12,
                                           face = "italic",
                                           rot = 270),
                         
  )
  
  
 
  ggsave(filename = name,plot = fig1,width = 7,height =7)
  #dev.off()
  
  
}


#make_figure(mod = results_sunny_bounded,name = "../bayesian_competition_ms//sunny_results.pdf")

#make_figure(mod = results_sunny_UNbounded,n_samples = 10,name = "../bayesian_competition_ms//sunny_results_unbounded.pdf")


make_figure(mod = results_woody_bounded,name = "../bayesian_competition_ms/woody_results.pdf")

#make_figure(mod = results_woody_UNbounded,n_samples = 10,name = "../bayesian_competition_ms/woody_results_unbounded.pdf")


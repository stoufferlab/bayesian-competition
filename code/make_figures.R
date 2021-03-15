library(ggplot2)
library(tidyverse)
source("code/gg_theme.R")
source("code/boxplot_proportions.R")


make_figure<-function(mod,
                      n_samples,
                      name){
  

  col1 <- rethinking::col.alpha("mediumseagreen", .8)
  col2 <- rethinking::col.alpha("grey50",.8)
  
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
      y = distance_growth_mean
    ), col= "#e3004e") +
    theme_alba +
    scale_color_manual(values = c(col2, col1)) +
    facet_grid(trcy_model~vero_model)  +
    geom_abline(intercept = 0, slope = 0, linetype ="dashed", col="grey50")
  #+
  # ylim(-1,1)+
  #xlim(0,0.05)
  
  fig1<- annotate_figure(integration,
                         
                         bottom = text_grob( expression('Feasibility domain,'~Omega),
                                             size = 12),
                         left = text_grob(expression('Distance from the boundary,'~theta),
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
  
  
  mod_proportions <- extract_proportion_samples(mod = mod,
                                                n_samples = n_samples)
  
  
  
  plot_proportions <-ggplot(mod_proportions) +
    geom_boxplot(mapping = aes(x=model,
                               y=proportions),
                 fill="mediumseagreen") +
    theme_alba +
    scale_x_discrete(labels=c("BH-BH","BH-LV","BH-RC","LV-BH","LV-LV","LV-RC","RC-BH","RC-LV","RC-RC")) 
  
  
  
  fig2 <- annotate_figure(plot_proportions,
                          bottom = text_grob("Model combination",
                                             size = 12),
                          left = text_grob("Proportion of coexistence",
                                           rot = 90, size =12))
  
  
  
  
  #pdf("test2", width = 7, height = 9)
  p<-ggarrange(fig1, fig2,
               ncol = 1,
               nrow = 2, heights = c(1,0.5))
  p
  
  ggsave(filename = name,plot = p,width = 7,height = 9)
  #dev.off()
  
  
}


make_figure(mod = example,n_samples = 100,name = "test.pdf")

library(ggplot2)
library(ggpubr)

# figures for the paper ---------------------------------------------------



col1 <- rethinking::col.alpha("mediumseagreen",0.3)
col2 <- rethinking::col.alpha("grey60",0.3)

col3<- rethinking::col.alpha("#e3004e", 1)

 
 
 
 open<- ggplot(results_sunny_bounded) +
   geom_point(
     mapping = aes(
       x = proportion,
       y = distance,
       col = as.factor(feasibility)
     ),
     show.legend = FALSE,
     size = 2
   ) +
   geom_point( mapping = aes(
     x = proportion_mean,
     y = distance_mean,
     shape= as.factor(feasibility_mean)),
     size=3,
     col=col3,
     show.legend = FALSE ) +
   scale_color_manual(values = c(col2, col1)) +
   facet_grid(trcy_model~vero_model)  +
   geom_abline(intercept = 0, slope = 0, linetype ="dashed", col="grey50")+
   xlim(0,1)+
   ylim(-25,25)+
   theme_bw()+
   theme(strip.background =element_rect(fill="white"))
 
 
 
 
 
 closed <- ggplot(results_woody_bounded) +
   geom_point(
     mapping = aes(
       x = proportion,
       y = distance,
       col = as.factor(feasibility)
     ),
     show.legend = FALSE,
     size = 2
   ) +
   geom_point( mapping = aes(
     x = proportion_mean,
     y = distance_mean,
     shape= as.factor(feasibility_mean)),
     size=3,
     col=col3,
     show.legend = FALSE ) +
   scale_color_manual(values = c(col2, col1)) +
   facet_grid(trcy_model~vero_model)  +
   geom_abline(intercept = 0, slope = 0, linetype ="dashed", col="grey50")+
   xlim(0,1)+
   ylim(-25,25)+
   theme_bw()+
   theme(strip.background =element_rect(fill="white"))
 
 
ggarrange(open, closed)
 
 
 
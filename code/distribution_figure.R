library(ggplot2)
library(ggpubr)

# figures for the paper ---------------------------------------------------

source("code/access_results.R")

sunny_results <- filter(final_results, environment=="open") %>% filter( vero_model != "Lotka-Volterra")%>% filter( trcy_model != "Lotka-Volterra" )
woody_results <- filter(final_results, environment=="woody")  %>% filter( vero_model != "Lotka-Volterra")%>% filter( trcy_model != "Lotka-Volterra" )

write_csv(sunny_results, file = "sunny.csv")
write_csv(woody_results, file = "woody.csv")


col1 <- rethinking::col.alpha("mediumseagreen",0.3)
col2 <- rethinking::col.alpha("grey60",0.3)
col3<- rethinking::col.alpha("#e3004e", 1)

xx <- expression("Relative coexistence ratio,"~rho)
yy <- expression("Distance from the edge,"~delta)


all_results<- rbind(sunny_results, woody_results) 

all_results <- all_results %>% unite("combo", vero_model:trcy_model, remove = FALSE)

 open<-ggplot(sunny_results) +
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
  xlim(0,2.5)+
  ylim(-20,10)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))+
  xlab(xx)+
  ylab(yy)

 
open <- annotate_figure(open,
                        top = text_grob( "Vellia rosea", size = 12, face = "italic"),
                        right = text_grob( "Trachymene cyanopetala", size = 12, face = "italic", rot=-90)
                        )
 
 
 ggsave(open, filename = "../bayesian_competition_ms/sunny_results.pdf", width = 8, height = 8)


# woody -------------------------------------------------------------------

 
 woody<-ggplot(woody_results) +
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
   theme_bw()+
   theme(strip.background =element_rect(fill="white"),
         axis.title = element_text(size = 12),
         strip.text = element_text(size = 12))+
   xlab(xx)+
   ylab(yy)+
   ylim(-20,10)+
   xlim(0,2.5)
 
 
 woody <- annotate_figure(woody,
                         top = text_grob( "Vellia rosea", size = 12, face = "italic"),
                         right = text_grob( "Trachymene cyanopetala", size = 12, face = "italic", rot=-90)
 )
 
 ggsave(woody, filename = "../bayesian_competition_ms/woody_results.pdf", width = 8, height = 8)
 

 

# push --------------------------------------------------------------------

all_results<- rbind(sunny_results, woody_results) 
 
 all_results <- all_results %>% unite("combo", vero_model:trcy_model, remove = FALSE)
 
 ggplot(all_results) +
   geom_point(aes(x=combo,
                  y=distance_mean,
                  col=environment),
              size=3)+
   theme(axis.text.x=element_text(angle=90,hjust=1))
 
 
 ggplot(all_results) +
   geom_point(aes(x=proportion_mean,
                  y=distance_mean,
                  col=as.factor(feasibility_mean),
                  shape=as.factor(combo)),
              size=4)+
   facet_grid(~environment)+
   scale_shape_manual(values = 0:9)+
   theme_bw()
 
 
 
 
 
 
 

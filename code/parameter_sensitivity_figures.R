library(ggplot2)
library(ggpubr)
library(tidyverse)


bh_bh_sensitiviyt <- read.csv("bh_bh_sensitiviyt.csv")

xx <- expression("C) Relative coexistence ratio,"~rho)
yy <- expression("D) Distance from the edge,"~delta)

parameters <- c( expression(alpha[ii]),
                 expression(alpha[ij]),
                 expression(alpha[ji]),
                 expression(alpha[jj]),
                 expression(lambda[i]),
                 expression(lambda[j])
                 )




sunny <- read.csv("sunny.csv")

sunny <- sunny %>% filter(
  vero_model == "Beverton-Holt",
  trcy_model == "Beverton-Holt"
)

woody <- read.csv("woody.csv")

woody <- woody %>% filter(
  vero_model == "Beverton-Holt",
  trcy_model == "Beverton-Holt"
)

results <- read.csv("results.csv")

results <- results %>% filter(
  vero_model == "Beverton-Holt",
  trcy_model == "Beverton-Holt"
)



area_feasible_mean <-
  data.frame(
    "area_feasible" = c(
      unique(sunny$area_mean),
      unique(woody$area_mean)
    ),
    "environment" = c("open", "woody")
  )

N <- nrow(bh_bh_sensitiviyt)/2

area_feasible_mean <- do.call("rbind", replicate(N, area_feasible_mean, simplify = FALSE))



area_alone_mean <-
  data.frame(
    "area_alone" = c(
      unique(sunny$area_alone_mean),
      unique(woody$area_alone_mean)
    ),
    "environment" = c("open", "woody")
  )


area_alone_mean <- do.call("rbind", replicate(N, area_alone_mean, simplify = FALSE))





proportion_mean <-
  data.frame(
    "proportion_mean_value" = c(
      unique(sunny$proportion_mean),
      unique(woody$proportion_mean)
    ),
    "environment" = c("open", "woody")
  )


proportion_mean <- do.call("rbind", replicate(N, proportion_mean, simplify = FALSE))






distance_mean <-
  data.frame(
    "distance_mean_value" = c(
      unique(sunny$distance_mean),
      unique(woody$distance_mean)
    ),
    "environment" = c("open", "woody")
  )

distance_mean <- do.call("rbind", replicate(N, distance_mean, simplify = FALSE))



# area_feasible -----------------------------------------------------------

area_together <-  ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=area_feasible),
               outlier.alpha = 0.5,
               outlier.size = 0.5)+
  theme_bw()+
  geom_abline(data = area_feasible_mean,
              aes(intercept = area_feasible, slope=0)
             )+
  facet_grid(~environment)+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size=12),
    axis.title = element_blank(),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  scale_x_discrete(labels=parameters)+
  ggtitle(  expression("A) Biologically-constrained feasibility domain,"~beta))
 




# area alone ----------------------------------------------------------------


area_alone <- ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x= parameter,y=area_alone),
               outlier.alpha = 0.5,
               outlier.size = 0.5)+
  geom_abline(data = area_alone_mean,
              aes(intercept = area_alone, slope=0)
  )+
  facet_grid(~environment)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size=12),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  scale_x_discrete(labels=parameters)+
  ggtitle( expression("B) Area in monocolture,"~gamma))




# proportion --------------------------------------------------------------


rho <- ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=proportion),
               outlier.alpha = 0.5,
               outlier.size = 0.5)+
  geom_abline(data = proportion_mean,
              aes(intercept = proportion_mean_value, slope=0)
  )+
  facet_grid(~environment)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  scale_x_discrete(labels=parameters)+
  ggtitle(xx)+
  ylim(0,2.5)

# distance ----------------------------------------------------------------


delta <-
  
  ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=distance),
               outlier.alpha = 0.5,
               outlier.size = 0.5
               )+
  geom_abline(data = distance_mean,
              aes(intercept = distance_mean_value, slope=0)
  )+
  facet_grid(~environment)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size=12),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  scale_x_discrete(labels=parameters)+
  ggtitle(yy)


# -------------------------------------------------------------------------

rho_delta <- ggarrange(
          rho, delta,
          ncol = 1,
          nrow = 2)



rho_delta<- annotate_figure(rho_delta,
                         bottom = text_grob( "Parameter", size = 15),
                         left  = text_grob( "Value", size = 15,rot=90)
)





areas <- ggarrange(
 area_together, area_alone,
  ncol = 1,
  nrow = 2)



areas<- annotate_figure(areas,
                            bottom = text_grob( "Parameter", size = 15),
                            left  = text_grob( "Value", size = 15,rot=90)
                        
)



together <- ggarrange(
  area_together, area_alone,
  rho, delta,
  align ="hv",
  ncol = 1,
  nrow = 4)



together <-  annotate_figure(together,
                             bottom = text_grob( "Parameter", size = 15),
                             left  = text_grob( "Value", size = 15,rot=90),
                           
                             
)


ggsave(together, filename = "../bayesian_competition_ms/params_results.pdf", width = 6, height = 10)

ggsave(areas, filename = "../bayesian_competition_ms/areas_results.pdf", width = 7, height = 6)




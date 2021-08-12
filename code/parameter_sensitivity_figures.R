library(ggplot2)
library(ggpubr)

bh_bh_sensitiviyt <- read_csv("bh_bh_sensitiviyt.csv")

xx <- expression("Relative coexistence ratio,"~rho)
yy <- expression("Distance from the edge,"~delta)

parameters <- c( expression(alpha[ii]),
                 expression(alpha[ij]),
                 expression(alpha[ji]),
                 expression(alpha[jj]),
                 expression(lambda[i]),
                 expression(lambda[j])
                 )




sunny <- read_csv("sunny.csv")

sunny <- sunny %>% filter(
  vero_model == "Beverton-Holt",
  trcy_model == "Beverton-Holt"
)

mean(sunny$area_alone)


# area_feasible -----------------------------------------------------------



area_together <-ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=area_feasible),
               fill)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size=15),
    axis.title = element_blank(),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  geom_abline(intercept = unique(sunny$area_mean),
              slope = 0,
              col = "grey50",
              linetype = "dashed",
              size = 1.3)+
  scale_x_discrete(labels=parameters)+
  ggtitle(  expression("Biologically-constrained feasibility domain,"~beta))




# area alone ----------------------------------------------------------------


area_alone <-ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=area_alone),
               fill)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size=15),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 15),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  geom_abline(intercept = unique(sunny$area_alone_mean),
              slope = 0,
              col = "grey50",
              linetype = "dashed",
              size = 1.3)+
  scale_x_discrete(labels=parameters)+
  ggtitle( expression("Area in monocolture,"~gamma))


# proportion --------------------------------------------------------------


rho <- ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=proportion),
               fill)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 15),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 15),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  geom_abline(intercept = unique(sunny$proportion_mean),
              slope = 0,
              col = "grey50",
              linetype = "dashed",
              size = 1.3)+
  scale_x_discrete(labels=parameters)+
  ggtitle(xx)

# distance ----------------------------------------------------------------


delta <-ggplot(bh_bh_sensitiviyt)+
  geom_boxplot(aes(x=parameter,
                   y=distance)
               )+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size=15),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 15),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  geom_abline(intercept = unique(sunny$distance_mean),
              slope = 0,
              col = "grey50",
              linetype = "dashed",
              size = 1.3)+
  scale_x_discrete(labels=parameters)+
  ggtitle(yy)


# -------------------------------------------------------------------------

ggarrange(area_together, area_alone,
          rho, delta,
          ncol = 2,
          nrow = 2)


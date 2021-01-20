library(ggplot2)
library(tidyverse)
source("code/gg_theme.R")

mod<- readRDS("~/bayesian-competition/Beverton-HoltandBeverton-Holt.RDS")


col1 <- rethinking::col.alpha("mediumseagreen", .3)
col2 <- rethinking::col.alpha("grey50",1)

saavedras <- ggplot(mod) +
  geom_point(
    mapping = aes(
      x = Omega_saaveda,
      y = theta_saavedra,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean_saaveda,
    y = theta_mean_saavedra
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  xlim(0,0.2)+
  ylim(0, 20)+
  ggtitle("Saavedra's Integration")



distance_2 <-ifelse(mod$feasibility ==1,mod$distance,mod$distance*-1)
mod$distance_2 <-distance_2

integration <- ggplot(mod) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = distance_2,
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
  xlim(0,0.2)+
  ylim(-4,4)+
  ggtitle("Bounded integration")

pdf("compare_integration.pdf", width = 10, height = 7)
ggarrange(saavedras,
          integration,
          ncol=2,
          nrow=1)
dev.off()
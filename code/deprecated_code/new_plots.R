source("code/gg_theme.R")
library(ggplot2)
library(tidyverse) 


results_sunny <- readRDS("results_sunny.RDS")

col1 <- rethinking::col.alpha("mediumseagreen", .5)
col2 <- rethinking::col.alpha("grey50", .5)




feasibility_plot_sunny<-ggplot(results_sunny) +
  geom_point(
    mapping = aes(
      x = Omega,
      y = theta_saavedra,
      col = as.factor(feasibility)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = theta_mean_saavedra
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  geom_abline(slope = 180, intercept = 0, col="grey50")+
  facet_grid(trcy_model~vero_model)

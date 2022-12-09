

library(ggplot2)

mod <- results_sunny_bounded 
mod2 <- results_woody_bounded 



col1 <- rethinking::col.alpha("mediumseagreen", .5)
col2 <- rethinking::col.alpha("grey60",.5)
col3<- rethinking::col.alpha("#e3004e", 1)



integration <- ggplot(mod) +
  geom_point(
    mapping = aes(
      x = Omega_saavedra,
      y = theta_saavedra,
      col = as.factor(feasibility_saavedra)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_saavedraa_mean,
    y = theta_saavedra_mean,
    shape= as.factor(feasibility_saavedra_mean)),
    size=3,
    col=col3,
    show.legend = FALSE ) +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_grid(trcy_model~vero_model) +
  xlim(0,0.25)




fig1<- annotate_figure(integration,
                       
                       bottom = text_grob( expression(Omega),
                                           size = 12),
                       
                       left = text_grob(expression(theta),
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


ggsave(filename = "../bayesian_competition_ms/saavedra_sunny.pdf",plot = fig1,width = 7,height =7)




# woody -------------------------------------------------------------------


integration <- ggplot(mod2) +
  geom_point(
    mapping = aes(
      x = Omega_saavedra,
      y = theta_saavedra,
      col = as.factor(feasibility_saavedra)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_saavedraa_mean,
    y = theta_saavedra_mean,
    shape= as.factor(feasibility_saavedra_mean)),
    size=3,
    col=col3,
    show.legend = FALSE ) +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_grid(trcy_model~vero_model) 




fig1<- annotate_figure(integration,
                       
                       bottom = text_grob( expression(Omega),
                                           size = 12),
                       
                       left = text_grob(expression(theta),
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


ggsave(filename = "../bayesian_competition_ms/saavedra_woody.pdf",plot = fig1,width = 7,height =7)









library(brms)
library(ggplot2)
library(bayesplot)

#We source everything known to human kind...

source("code/gg_theme.R")
source("code/read_models.R")
source("code/model_toolbox.R")


#isaacs mean estimates
gi<- 0.9641188
si<- 0.9654804
gj<- 0.4743974
sj<- 0.9693324


multiple_models <-function(..., s,g) {
  
  models <- list(...)
  posteriors <- lapply(models, function(m) {
    one_posterior <- posterior_parameters(m,s,g)
    return(one_posterior)
  })
  all_posteriors <- do.call(rbind, posteriors)
  
}


posterior_distributions <- function(model){

  model_posterior <- multiple_models(model,
                                     s=1,
                                     g=1)
  
  name <- unique(model_posterior$model)
  

  competition  <- model_posterior %>% select(
                               alphaii,
                               alphaii_env,
                               alphaij,
                               alphaij_env,
                               alphaik,
                               alphaik_env) 

  colnames(competition) <- c("Intraspecific",
                             "Intraspecific woody",
                             "Interespecific",
                             "Interespecific woody",
                             "Other",
                             "Other woody")
  
  
competition_plot<-  mcmc_areas(competition,
             prob = 0.8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(size = 12),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  xlim(-0.03, 0.16)
  

seeds  <- model_posterior %>% select(lambda,
                                     lambda_env)

colnames(seeds) <- c("Fecundity",
                           "Fecundity woody")


seed_plot<-  mcmc_areas(seeds,
                               prob = 0.8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(size = 12),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  xlim(0, 13)

together <- ggarrange(competition_plot,
          seed_plot,
          nrow = 2, 
          heights = c(1,0.5))

together <- annotate_figure(together, top = text_grob(name, size = 15)
                            )  

return(together)
}



  

  
# vero --------------------------------------------------------------------


vero_bh <- posterior_distributions(model = vero_bh_multispecies_poisson.rds)


vero_rc <- posterior_distributions(model = vero_rc_multispecies_poisson.rds)



vero_dist <- ggarrange(vero_bh,
                       vero_rc,
                       ncol = 2
                       )


vero_dist <- annotate_figure(
  vero_dist,
  top = text_grob(
    "Vellia rosea",
    face = "italic",
    size = 15
  ),
  bottom = text_grob("Value", size = 15),
  left = text_grob("Parameter", size = 15, rot = 90),
  fig.lab.pos = "top.right"
)


pdf(file = "~/bayesian_competition_ms/vero_dist.pdf", height = 8, width = 8)

vero_dist

dev.off()




# trcy --------------------------------------------------------------------




posterior_distributions <- function(model){
  
  model_posterior <- multiple_models(model,
                                     s=1,
                                     g=1)
  
  name <- unique(model_posterior$model)
  
  
  competition  <- model_posterior %>% select(
    alphaii,
    alphaii_env,
    alphaij,
    alphaij_env,
    alphaik,
    alphaik_env) 
  
  colnames(competition) <- c("Intraspecific",
                             "Intraspecific woody",
                             "Interespecific",
                             "Interespecific woody",
                             "Other",
                             "Other woody")
  
  
  competition_plot<-  mcmc_areas(competition,
                                 prob = 0.8) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      strip.text = element_text(size = 12),
      legend.text.align = 0,
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlim(-0.03, 0.23)
  
  
  seeds  <- model_posterior %>% select(lambda,
                                       lambda_env)
  
  colnames(seeds) <- c("Fecundity",
                       "Fecundity woody")
  
  
  seed_plot<-  mcmc_areas(seeds,
                          prob = 0.8) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      strip.text = element_text(size = 12),
      legend.text.align = 0,
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlim(0, 43)
  
  together <- ggarrange(competition_plot,
                        seed_plot,
                        nrow = 2, 
                        heights = c(1,0.5))
  
  together <- annotate_figure(together, top = text_grob(name, size = 15)
  )  
  
  return(together)
}



trcy_bh <- posterior_distributions(model = trcy_bh_multispecies_poisson.rds)

trcy_rc <- posterior_distributions(model = trcy_rc_multispecies_poisson.rds)



trcy_dist <- ggarrange(trcy_bh,
                       trcy_rc,
                       ncol = 2
)

trcy_dist <- annotate_figure(
  trcy_dist,
  top = text_grob(
    "Trachymene cyanopetala",
    face = "italic",
    size = 15),
  bottom = text_grob("Value", size = 15),
  left = text_grob("Parameter", size = 15, rot = 90),
  fig.lab.pos = "top.right"
)

pdf(file = "~/bayesian_competition_ms/trcy_dist.pdf", height = 8, width = 8)

trcy_dist

dev.off()






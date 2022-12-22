
library(bayesplot)
library(ggplot2)
require(ggpubr)
require(tidyverse)
require(here)

source(here("code/lib/load_models.R"))
source(here("code/lib/posterior_utils.R"))

# function that for a model generates plots of the posteriors of key derived parameters
posterior_distributions <- function(model, xlim1=NULL, xlim2=NULL){

  model_posterior <- posterior_parameters(model,s=NA,g=NA,model$name)

  if(model$name == "Null"){
    # treat the posterior as a bunch of zeros
    competition  <- data.frame(
        alphaii_open  = rep(0, 10),
        alphaii_woody = rep(0, 10),
        alphaij_open  = rep(0, 10),
        alphaij_woody = rep(0, 10),
        alphaik_open  = rep(0, 10),
        alphaik_woody = rep(0, 10)
      )
  }else{
    competition  <- model_posterior %>%
      select(
        alphaii_open,
        alphaii_woody,
        alphaij_open,
        alphaij_woody,
        alphaik_open,
        alphaik_woody
      )
  }

  colnames(competition) <- c(" Intraspecific (open)",
                             "Intraspecific (woody)",
                             " Interspecific (open)",
                             "Interspecific (woody)",
                             "         Other (open)",
                             "        Other (woody)")
  
  competition_plot <- mcmc_intervals(
      competition,
      prob = 0.5
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis.title = element_blank(),
      strip.text = element_text(size = 12),
      legend.text.align = 0,
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab('Value') +
    theme(axis.title.x = element_text(size=15))

  if(!is.null(xlim1)){
    competition_plot <- competition_plot + xlim(xlim1[1], xlim1[2])
  }

  seeds  <- model_posterior %>% select(lambda_open,
                                       lambda_woody)

  colnames(seeds) <- c(
    "     Fecundity (open)",
    "    Fecundity (woody)"
  )

  seed_plot <- mcmc_intervals(
      seeds,
      prob = 0.5
    ) +
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
    )

  if(!is.null(xlim2))
    seed_plot <- seed_plot + xlim(xlim2[1], xlim2[2])

  together <- ggarrange(
    seed_plot,
    competition_plot,
    nrow = 2, 
    heights = c(0.25,0.75)
  )

  together <- annotate_figure(
    together,
    top = text_grob(model$name, size = 15)
  )  

  return(together)
}
  
# vero --------------------------------------------------------------------

vero_nl <- posterior_distributions(
  model = focal_models[["vero"]][["interaction_free"]],
  xlim1 = c(-0.1, 0.1),
  xlim2 = c(0, 20)
)

# vero_lv <- posterior_distributions(
#   model = focal_models[["vero"]][["lotka_volterra"]],
#   # xlim1 = c(-0.1, 0.4),
  # xlim2 = c(0, 30)
# )

vero_bh <- posterior_distributions(
  model = focal_models[["vero"]][["beverton_holt"]],
  xlim1 = c(-0.2, 0.2),
  xlim2 = c(0, 20)
)

vero_rc <- posterior_distributions(
  model = focal_models[["vero"]][["ricker"]],
  xlim1 = c(-0.1, 0.1),
  xlim2 = c(0, 20)
)

vero_dist <- ggarrange(
  vero_nl,
  # vero_lv,
  vero_rc,
  vero_bh,
  ncol = 3
)

vero_dist <- annotate_figure(
  vero_dist,
  top = text_grob(
    "Goodenia rosea",
    face = "italic",
    size = 15
  ),
  # bottom = text_grob("Value", size = 15),
  left = text_grob("Parameter", size = 15, rot = 90),
  fig.lab.pos = "top.right"
)

pdf(file = here("figures/goro_dist.pdf"), height = 8, width = 14)
vero_dist
dev.off()

# trcy --------------------------------------------------------------------

trcy_nl <- posterior_distributions(
  model = focal_models[["trcy"]][["interaction_free"]],
  xlim1 = c(-0.1, 0.1),
  xlim2 = c(0, 80)
)

# trcy_lv <- posterior_distributions(
#   model = focal_models[["trcy"]][["lotka_volterra"]],
#   # xlim1 = c(-0.15, 0.3),
#   # xlim2 = c(0, 80)
# )

trcy_bh <- posterior_distributions(
  model = focal_models[["trcy"]][["beverton_holt"]],
  xlim1 = c(-0.5, 0.5),
  xlim2 = c(0, 80)
)

trcy_rc <- posterior_distributions(
  model = focal_models[["trcy"]][["ricker"]],
  xlim1 = c(-0.08, 0.08),
  xlim2 = c(0, 80)
)

trcy_dist <- ggarrange(
  trcy_nl,
  # trcy_lv,
  trcy_rc,
  trcy_bh,
  ncol = 3
)

trcy_dist <- annotate_figure(
  trcy_dist,
  top = text_grob(
    "Trachymene cyanopetala",
    face = "italic",
    size = 15),
  # bottom = text_grob("Value", size = 15),
  left = text_grob("Parameter", size = 15, rot = 90),
  fig.lab.pos = "top.right"
)

pdf(file = here("figures/trcy_dist.pdf"), height = 8, width = 14)
trcy_dist
dev.off()

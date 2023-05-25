
# define functions used below
require(bayesplot)
require(brms)
require(here)
require(ggplot2)
require(ggpubr)

# general code that finds and loads all candidate model fits
source(here("code/lib/load_models.R"))

# to make sure we always get the same pp_check
set.seed(123456)

bh <- pp_check(focal_models[["vero"]][["beverton_holt"]], discrete=TRUE, type="ecdf_overlay") +
        labs(x="Number of seeds", y="Cumulative distribution")
rick <- pp_check(focal_models[["vero"]][["ricker"]], discrete=TRUE, type="ecdf_overlay") +
        labs(x="Number of seeds", y="Cumulative distribution")
null <- pp_check(focal_models[["vero"]][["interaction_free"]], discrete=TRUE, type="ecdf_overlay") +
        labs(x="Number of seeds", y="Cumulative distribution")

# output filename
filename <- paste0(
    "figures/posterior_predictive_checks_goro.pdf"
)

together <- bayesplot::bayesplot_grid(
    null,
    rick,
    bh,
    grid_args=list(ncol=3),
    legends=FALSE,
    subtitles=c("Null","Ricker","Beverton-Holt")
)

together <- annotate_figure(
  together,
  top = text_grob(
    "Goodenia rosea",
    face = "italic",
    size = 15
  ),
  # bottom = text_grob("Value", size = 15),
  # left = text_grob("Parameter", size = 15, rot = 90),
  fig.lab.pos = "top.right"
)

ggsave(filename = here::here(filename), plot=together, height=3)

bh <- pp_check(focal_models[["trcy"]][["beverton_holt"]], discrete=TRUE, type="ecdf_overlay") +
        labs(x="Number of seeds", y="Cumulative distribution")
rick <- pp_check(focal_models[["trcy"]][["ricker"]], discrete=TRUE, type="ecdf_overlay") +
        labs(x="Number of seeds", y="Cumulative distribution")
null <- pp_check(focal_models[["trcy"]][["interaction_free"]], discrete=TRUE, type="ecdf_overlay") +
        labs(x="Number of seeds", y="Cumulative distribution")

# output filename
filename <- paste0(
    "figures/posterior_predictive_checks_trcy.pdf"
)

together <- bayesplot::bayesplot_grid(
    null,
    rick,
    bh,
    grid_args=list(ncol=3),
    legends=FALSE,
    subtitles=c("Null","Ricker","Beverton-Holt")
)

together <- annotate_figure(
  together,
  top = text_grob(
    "Trachymene cyanopetala",
    face = "italic",
    size = 15
  ),
  # bottom = text_grob("Value", size = 15),
  # left = text_grob("Parameter", size = 15, rot = 90),
  fig.lab.pos = "top.right"
)

ggsave(filename = here::here(filename), plot=together, height=3)

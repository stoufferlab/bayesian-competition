
library(tidyverse)
library(xtable)
require(here)

source(here("code/lib/load_models.R"))
source(here("code/lib/compare_models.R"))

# vero

# vero model comparison in open environment / woody==0
vero_open <- model_compare(
  focal_models[["vero"]][["interaction_free"]],
  focal_models[["vero"]][["beverton_holt"]],
  focal_models[["vero"]][["ricker"]],
  woody=0
)
colnames(vero_open) <- paste0(colnames(vero_open),"_open")

# vero model comparison in woody environment / woody==1
vero_woody <- model_compare(
  focal_models[["vero"]][["interaction_free"]],
  focal_models[["vero"]][["beverton_holt"]],
  focal_models[["vero"]][["ricker"]],
  woody=1
)
colnames(vero_woody) <- paste0(colnames(vero_woody),"_woody")

# vero model comparison in both environments
vero_all <- model_compare(
  focal_models[["vero"]][["interaction_free"]],
  focal_models[["vero"]][["beverton_holt"]],
  focal_models[["vero"]][["ricker"]],
  woody=NULL
)
colnames(vero_all) <- paste0(colnames(vero_all),"_all")

vero <- vero_open %>%
  rownames_to_column(var="Model") %>%
  left_join(
    rownames_to_column(vero_woody, var="Model"),
    by="Model"
  ) %>%
  left_join(
    rownames_to_column(vero_all, var="Model"),
    by="Model"
  )
vero$Species <- "Goodenia_rosea"

# trcy

# trcy model comparison in open environment / woody==0
trcy_open <- model_compare(
  focal_models[["trcy"]][["interaction_free"]],
  focal_models[["trcy"]][["beverton_holt"]],
  focal_models[["trcy"]][["ricker"]],
  woody=0
)
colnames(trcy_open) <- paste0(colnames(trcy_open),"_open")

# trcy model comparison in woody environment / woody==1
trcy_woody <- model_compare(
  focal_models[["trcy"]][["interaction_free"]],
  focal_models[["trcy"]][["beverton_holt"]],
  focal_models[["trcy"]][["ricker"]],
  woody=1
)
colnames(trcy_woody) <- paste0(colnames(trcy_woody),"_woody")

# trcy model comparison in both environments
trcy_all <- model_compare(
  focal_models[["trcy"]][["interaction_free"]],
  focal_models[["trcy"]][["beverton_holt"]],
  focal_models[["trcy"]][["ricker"]],
  woody=NULL
)
colnames(trcy_all) <- paste0(colnames(trcy_all),"_all")

trcy <- trcy_open %>%
  rownames_to_column(var="Model") %>%
  left_join(
    rownames_to_column(trcy_woody, var="Model"),
    by="Model"
  ) %>%
  left_join(
    rownames_to_column(trcy_all, var="Model"),
    by="Model"
  )
trcy$Species <- "Trachymene_cyanopetala"

# both focal species together

values <- rbind(vero,trcy)

# DEBUG: write this to a file?
xtable(values)

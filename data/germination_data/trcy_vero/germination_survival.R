rm(list = ls())
library(tidyverse)
library(lme4)

### for the fill estimates
data_fill <- read_csv("data/trcy_vero_fill.csv") %>%
  mutate(percentage.filled = ((no.seeds.visible - no.seeds.unfilled) / no.seeds.visible)) %>%
  select(plot,
         species,
         no.seeds.visible,
         no.seeds.unfilled,
         percentage.filled)

trcy_fill <- data_fill %>% filter(species == "Trachymene cyanopetala")
vero_fill <- data_fill %>% filter(species == "Velleia rosea")


trcy_fill_model<-glmer(cbind(no.seeds.visible-no.seeds.unfilled, (no.seeds.unfilled)) ~1 + (1|plot), 
                       family=binomial, 
                       data = trcy_fill)

vero_fill_model<-glmer(cbind(no.seeds.visible-no.seeds.unfilled, (no.seeds.unfilled)) ~1 + (1|plot), 
                       family=binomial, 
                       data = vero_fill)


trcy_fill_estimate <- fixef(trcy_fill_model) %>% plogis()
vero_fill_estimate <- fixef(vero_fill_model) %>% plogis()






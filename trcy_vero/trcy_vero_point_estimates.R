rm(list = ls())
library(tidyverse)
library(lme4)
####bring in germiantion data
trcy_vero<-read_csv("data/trcy_vero.csv") %>%
  dplyr::select(species, plot, light_treatment = `light/ dark`, plate_number = 'plate #',
                num_seeds = '# seeds', emerged_inviable = 'emerged-inviable',
                nongerm_filled = 'nongerminant-filled', nongerm_unfilled = 'nongerminant-unfilled',
                no_radical = 'cotyledons only',
                week_1 = 'week 1', week_2 = 'week 2', week_3 = 'week 3', week_4 = 'week 4', week_5 = 'week 5', cut_test_performed) %>%
  dplyr::mutate(prop_germ = week_5/num_seeds)

####seed fill data
trcy_vero_fill<-read_csv("data/trcy_vero_fill.csv")%>%
  mutate(percentage.filled=((no.seeds.visible-no.seeds.unfilled)/no.seeds.visible))%>%
  select(plot, species,no.seeds.visible,no.seeds.unfilled,percentage.filled)

#####calculate the point estiamte and CI for trcy seed fill first
trcy.fill<-glmer(cbind(no.seeds.visible-no.seeds.unfilled, (no.seeds.unfilled)) ~1 + (1|plot), family=binomial, trcy_vero_fill[trcy_vero_fill$species=="Trachymene cyanopetala",])

####plogis back-transform point estimate + and - standard error
trcy.mean.fill<-plogis(3.4534)
trcy.fill.lower.ci<-plogis(3.4534-0.2118)
trcy.fill.upper.ci<-plogis(3.4534+0.2118)

#####calculate the point estiamte and CI for vero seed fill next
vero.fill<-glmer(cbind(no.seeds.visible-no.seeds.unfilled, (no.seeds.unfilled)) ~1 + (1|plot), family=binomial, trcy_vero_fill[trcy_vero_fill$species=="Velleia rosea",])

####plogis back-transform point estimate + and - standard error
vero.mean.fill<-plogis(3.3311)
vero.fill.lower.ci<-plogis(3.3311-0.2343)
vero.fill.upper.ci<-plogis(3.3311+0.2343)


###join fill to germination
trcy_vero<-left_join(trcy_vero, select(trcy_vero_fill,plot, species, percentage.filled,no.seeds.visible,no.seeds.unfilled))

###adjust the number of seeds in petri dishes by plot-specific fill rates
trcy_vero$adjusted.num.seeds=round(trcy_vero$num_seeds*trcy_vero$percentage.filled, digits=0)

#####function to bump up the adjusted seeds if high seed fill in xray conflicted with the adjusted number of seeds in petri dish
trcy_vero$adjusted.num.seeds<-ifelse(with(trcy_vero, adjusted.num.seeds-no_radical-emerged_inviable-week_5)<0, with(trcy_vero, no_radical+emerged_inviable+week_5), trcy_vero$adjusted.num.seeds)

#####first just calculate the proportion of germination in each petri dish
trcy_vero$prop_germ<-trcy_vero$week_5/(trcy_vero$adjusted.num.seeds-trcy_vero$no_radical-trcy_vero$emerged_inviable)


#####just extract the treatment and plot where species germianted best (i.e light or dark for each plot)
higher.germ<-trcy_vero%>% 
  group_by(plot,species)%>% 
  slice(which.max(prop_germ))

#####calculate the point estiamte and CI for trcy gemriantion first
trcy.germ<-glmer(cbind(week_5, ((adjusted.num.seeds-no_radical-emerged_inviable)-week_5)) ~1 + (1|plate_number), family=binomial, higher.germ[higher.germ$species=="Trachymene cyanopetala",])

####plogis back-transform point estimate + and - standard error
trcy.mean.germ<-plogis(-0.1025)
trcy.germ.lower.ci<-plogis(-0.1025-0.2442)
trcy.germ.upper.ci<-plogis(-0.1025+0.2442)

#####calculate the point estiamte and CI for vero gemriantion second
vero.germ<-glmer(cbind(week_5, ((adjusted.num.seeds-no_radical-emerged_inviable)-week_5)) ~1 + (1|plate_number), family=binomial, higher.germ[higher.germ$species=="Velleia rosea",])

####plogis back-transform point estimate + and - standard error
vero.mean.germ<-plogis(3.2910)
vero.germ.lower.ci<-plogis(3.2910-0.4407)
vero.germ.upper.ci<-plogis(3.2910+0.4407)




# table -------------------------------------------------------------------


results <- data.frame("Species"= c("Vellia rosea", "Trachymene cyanopetala"),
                      "Germination"= c(vero.mean.germ, trcy.mean.germ) %>% round(digits = 3),
                      "Survival"=c(vero.mean.fill, trcy.mean.fill) %>% round(digits = 3))


tab <- xtable::xtable(results)





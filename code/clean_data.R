library(tidyverse)
#the plde and vero data set

dat2<-read.csv("data/VEROvTRCY.FinalDatasetFromAIP.csv")

#dat1<-as_tibble(dat1)
dat2<-as_tibble(dat2)


#rename in a way that is easy to understand
vero_trcy_data <- dat2 %>%  rename( verotreatment = Nv, veroplanted=Veroseeds.planted, verodensity=Nv.1, trcytreatment=Nt, trcyplanted=TRCYseeds.planted, trcydensity=Nt.1, originalenv=Original.Env, exposedenv=Exp.Env, totalseeds=number.of.seeds)

#make the  exposed environment a variable 
vero_trcy <- vero_trcy_data %>% mutate( env = ifelse(exposedenv=="woody",1,0 )) %>% drop_na()

#separate by species
vero_focal<- vero_trcy %>% filter(focal == "V") %>% rename("conspecifics" = verodensity, "heterospecifics" = trcydensity )
trcy_focal<- vero_trcy %>% filter(focal == "T") %>% rename("conspecifics" = trcydensity,  "heterospecifics"= verodensity)



rm(dat2,vero_trcy_data, vero_trcy)

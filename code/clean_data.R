library(tidyverse)
#the plde and vero data set
dat1<-read.csv("data/PLDEvVero.FinalDatasetFromAIP.csv")
#the trycy data set
dat2<-read.csv("data/VEROvTRCY.FinalDatasetFromAIP.csv")

dat1<-as_tibble(dat1)
dat2<-as_tibble(dat2)


#rename variables so it is easy to understand them
# vero_plde<- dat1 %>% rename(verotreatment=Vtrt, veroplanted=VEROseeds, verodensity=Vn, pldetreatment=Ptrt, pldeplanted=PLDEseeds, pldedensity=Pn, totalother=total.other, originalenv=Original.Env, exposedenvo=Exp.Env, focal=Focal ,totalseeds=Seeds.total)
# #make the exposed environment a variable
# vero_plde<- vero_plde %>% mutate( env = ifelse(originalenv=="shade",1,0))
# 
# 
# 

 
#rename in a way that is easy to understand
vero_trcy_data <- dat2 %>%  rename( verotreatment = Nv, veroplanted=Veroseeds.planted, verodensity=Nv.1, trcytreatment=Nt, trcyplanted=TRCYseeds.planted, trcydensity=Nt.1, originalenv=Original.Env, exposedenv=Exp.Env, totalseeds=number.of.seeds)

#make the  exposed environment a variable 
vero_trcy <- vero_trcy_data %>% mutate( env = ifelse(originalenv=="woody",1,0 )) %>% drop_na()




rm(dat1,dat2, vero_trcy_data)

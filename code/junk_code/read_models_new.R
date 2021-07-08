library(brms)
library(tidyverse)
#we go where the models are
setwd("./model_objects_new")
#we list them
files <-list.files(pattern = ".rds") 
#we read them
all_models<-lapply(files,readRDS)
#we name them
names(all_models)<-files
#we upload them to the environment
list2env(all_models , envir = .GlobalEnv)

#we remove what we do not want
rm(all_models, files)


#we add a name to each model and the constraints to its growth rate
vero_bh_multispecies_poisson.rds$name<-"Beverton-Holt"
vero_bh_multispecies_poisson.rds$constraints <- c(-1, Inf)

vero_lv_multispecies_poisson.rds$name<- "Lotka-Volterra"
vero_lv_multispecies_poisson.rds$constraints<- c(-Inf,1)

vero_rc_multispecies_poisson.rds$name<- "Ricker"
vero_rc_multispecies_poisson.rds$constraints<- c(-Inf, Inf)



trcy_bh_multispecies_poisson.rds$name<-"Beverton-Holt"
trcy_bh_multispecies_poisson.rds$constraints<-c(-1, Inf)

trcy_lv_multispecies_poisson.rds$name<- "Lotka-Volterra"
trcy_lv_multispecies_poisson.rds$constraints<- c(-Inf,1)


trcy_rc_multispecies_poisson.rds$name<- "Ricker"
trcy_rc_multispecies_poisson.rds$constraints<- c(-Inf, Inf)



#we leave
setwd("..")
print(getwd())
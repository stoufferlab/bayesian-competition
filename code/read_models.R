library(brms)
library(tidyverse)
#we go where the models are
setwd("./model_objects")
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

#we leave
setwd("..")
print(getwd())
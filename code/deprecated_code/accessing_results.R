library(tidyverse)

files <-list.files(pattern = ".RDS") 
#we read them
all_results<-lapply(files,readRDS)
#we name them
names(all_results)<-files

#we upload them to the environment
list2env(all_results , envir = .GlobalEnv)

models_names <- lapply(seq_along(all_results), function(i,all_results){
     
      full_name <-   names(all_results)[[i]]
        name <- strsplit(full_name, ".RDS")
        species_name<- strsplit( name[[1]], "and")
        vero_name <- species_name[[1]][1]
        trcy_name <- species_name[[1]][2]
         all_results[[i]][,"vero_model"] <-  vero_name
         all_results[[i]][,"trcy_model"] <-  trcy_name
       return(all_results[[i]])  


}, all_results = all_results)

all_together <- do.call(rbind,models_names)

saveRDS(all_together, file = "results_sunny.RDS")

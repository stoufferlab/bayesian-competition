library(brms)
library(tidyverse)
#we go where the models are
setwd("./integration_objects")
#we list them
files <-list.files(pattern = ".RDS")
#we read them
all_models<-lapply(files,readRDS)


results <- lapply(seq_along(all_models), function(i){
   name <- files[i]
   
   decomposed_name <-stringr::str_split(name, pattern = "_")[[1]][3]
   
   mod <- all_models[[i]]
  
  if(decomposed_name=="FALSE.RDS"){
    print("k")
    mod$environment <- "open"
  }else{
    mod$environment <- "woody"
  }
  return(mod)
  
})


final_results <- do.call(rbind,results)
#we leave
setwd("..")

write_csv(final_results, file = "results.csv")


print(getwd())
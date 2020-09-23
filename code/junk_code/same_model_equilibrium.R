library(brms)
library(ggplot2)

 source("code/growth_rates.R")
 source("code/model_equilibrium.R")
 source("code/model_toolbox.R")
gi<-.372
si<-.556
gj<-.258
sj<-.033




multiple_equilibriums<-function(models,model_names,growth_fucntions, equilibrium_functions, s,g){
 
  all_posteriors<-data.frame()
   
  
   for(i in 1:length(models)){
     if(model_names[i]=="Hassell"){ exp_param= TRUE}else{exp_param=FALSE}
    posterior <- posterior_parameters(model=models[[i]],growth_fun = growth_fucntions[[i]],equilibrium_fun = equilibrium_functions[[i]],s = si,g = gi,exp_param = exp_param)
    posterior$model <- model_names[i]
    all_posteriors<-rbind(all_posteriors,posterior)
  }
  

  ggplot(all_posteriors)+
    geom_density(mapping = aes(x = equilibrium, fill=model) ,alpha=0.8)+
    scale_fill_manual(values=c("#109877", "#109798"))+
    theme_bw()+
    xlab("Monoculture equilibrium")+
    xlim(0,1000)

# 
#   ggplot(all_posteriors)+
#     geom_density(mapping = aes(x = env_equilibrium, fill=model) ,alpha=0.8)+
#     scale_fill_manual(values=c("#109877", "#109798"))+
#     theme_bw()+
#     xlab("Monoculture equilibrium woody environment") +
#     xlim(0,1000)
#   
  
    #scale_fill_manual(values=c("#109877", "#109798", "#107598", "#105398", "111098"))
    
}

#order is important in this function, each model has its own growth function and equilibrium function, and name
vero_models <-
  list(vero_bh_multispecies_zippoisson,
       vero_rc_multispecies_zippoisson)
trcy_models <- list( trcy_bh_multispecies_zipoisson)

model_names <- c("Beverton-Holt", "Ricker")
growth_functions <- list(bh_growth, rc_growth, lv_growth, hs_growth)
equilibrium_functions <-
  list(bh_equilibrium,
       rc_equilibrium,
       lv_equilibrium,
       hs_equilibrium)

vero<-multiple_equilibriums( models = vero_models,
                       model_names = model_names,
                       growth_fucntions = growth_functions,
                       equilibrium_functions = equilibrium_functions,
                       s = si,
                       g = gi)


trcy<-multiple_equilibriums( models = trcy_models,
                             model_names = model_names,
                             growth_fucntions = growth_functions,
                             equilibrium_functions = equilibrium_functions,
                             s = si,
                             g = gi)




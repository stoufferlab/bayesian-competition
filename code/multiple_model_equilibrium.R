library(brms)
library(ggplot2)
library(ggpubr)

#We source everything known to human kind...
source("code/growth_rates.R")
source("code/gg_theme.R")
source("code/read_models.R")
source("code/model_toolbox.R")


#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033



#A function that returns a data frame with posterior equilibrium values ofr each models
#modesl is a list of models of the same species
#names is the names of each model, IN THE SAME ORDER as in the list models
#growth functions are the growth functions of each model, in the same order as in the list models
#s and g are the respective survival and germination of each species
multiple_equilibriums <-
  function(models,
           model_names,
           growth_fucntions,
           s,
           g) {
    all_posteriors <- data.frame()
    
    
    for (i in 1:length(models)) {
      if (model_names[i] == "Hassell") {
        exp_param = TRUE
      } else{
        exp_param = FALSE
      }
      posterior <-
        posterior_parameters(
          model = models[[i]],
          growth_fun = growth_fucntions[[i]],
          s = s,
          g = g,
          exp_param = exp_param
        )
      posterior$model <- model_names[i]
      all_posteriors <- rbind(all_posteriors, posterior)
    }
    
  return(all_posteriors)
    
  }

#order is important in this function, each model has its own growth function and equilibrium function, and name


#This is the order : beverton holt (bh), lotka volterra (lv), ricker (rc) and hassell (hs)
model_names <- c("Beverton-Holt", "Lotka-Volterra", "Ricker", "Hassell")
growth_functions <- list(bh_growth, lv_growth, rc_growth, hs_growth)


vero_models <-list(vero_bh_multispecies_poisson.rds,
                   vero_lv_multispecies_poisson.rds,
                   vero_rc_multispecies_poisson.rds,
                   vero_hs_multispecies_poisson.rds)

trcy_models<-list(trcy_bh_multispecies_poisson.rds,
                  trcy_lv_multispecies_poisson.rds,
                  trcy_rc_multispecies_poisson.rds,
                  trcy_hs_multispecies_poisson.rds)



vero <- multiple_equilibriums(
  models = vero_models,
  model_names = model_names,
  growth_fucntions = growth_functions,
  s = si,
  g = gi
)


trcy <- multiple_equilibriums(
  models = trcy_models,
  model_names = model_names,
  growth_fucntions = growth_functions,
  s = sj,
  g = gj
)

trcy$model<-factor(trcy$model, levels = model_names)
vero$model<-factor(vero$model, levels= model_names)


  
vero_0<-ggplot(vero) +
  geom_density(mapping = aes(x = equilibrium, fill = model, linetype=model) ,
               alpha = 0.8, show.legend = TRUE) +
  scale_linetype_manual(values=c("solid","dashed","twodash", "dotted"))+
  scale_fill_manual(values=palette_alba)+
  theme_alba+
  xlim(0,400) +
  ylim(0,0.027)



trcy_0<-ggplot(trcy) +
  geom_density(mapping = aes(x = equilibrium, fill = model, linetype=model) ,
               alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values=palette_alba)+
  theme_alba +
  xlim(0,400)+
  ylim(0,0.027)



vero_1<-ggplot(vero) +
  geom_density(mapping = aes(x = env_equilibrium, fill = model, linetype=model) ,
               alpha = 0.8,show.legend = TRUE) +
  scale_fill_manual(values=palette_alba) +
  theme_alba +
  xlim(0,1500)+
  ylim(0,0.0125)


trcy_1<-ggplot(trcy) +
  geom_density(mapping = aes(x = env_equilibrium, fill = model, linetype=model) ,
               alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values=palette_alba)+
  theme_alba + 
  xlim(0,1500) +
  ylim(0,0.0125)







all_figures<-ggarrange( vero_0, trcy_0 ,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)


env_figures<-ggarrange( vero_1, trcy_1 ,
                        labels = c("A", "B"),
                        ncol = 1, nrow = 2)


pdf(file = "results/equilibrium.pdf", width = 6, height = 7.5/1.6)
annotate_figure(all_figures,
             
                bottom = text_grob("Monoculture equilibrium abundance",
                                  size = 14),
                left = text_grob("Density",  rot = 90, size = 14),
                
)
 
dev.off()





pdf(file = "results/env_equilibrium.pdf", width = 6, height = 7.5/1.6)
annotate_figure(env_figures,
                
                bottom = text_grob("Monoculture equilibrium abundance",
                                   size = 14),
                left = text_grob("Density",  rot = 90, size = 14),
                
)

dev.off()
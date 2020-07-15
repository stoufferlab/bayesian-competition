#test for figure regarding theta vs. omega


source("code/model_toolbox.R")

fixed_coexistence(BEV_vero,BEV_trcy,"Beverton-Holt ")
pp<-posterior_feasibility(BEV_vero,BEV_trcy) 

pp$feasibility_results<-as.factor(pp$feasibility_results)


col1<-col.alpha("deepskyblue4",1)
col2<-col.alpha("firebrick1",.3)
ggplot(pp) + geom_point(mapping = aes(x=omega_results,y=theta_results, col= feasibility_results)) +  
  scale_color_manual(values = c(col2,col1)) +
  labs(col = "Coexistence") + 
  xlab("Omega (niche differences)")+
  ylab("Theta (fitness differences)")+
  theme_classic()
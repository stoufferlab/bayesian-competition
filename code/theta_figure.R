#test for figure regarding theta vs. omega
require(gridExtra)
source("code/model_toolbox.R")

make_figure_theta <- function(vero_model,vero_function, vero_exp,trcy_model,trcy_function,trcy_exp,env){
 
   col1<-rethinking::col.alpha("dodgerblue3",.3)
   col2<-rethinking::col.alpha("firebrick3",.3)
  
  post        <-  feasibility_wrapper(vero_model,vero_function, vero_exp,trcy_model,trcy_function,trcy_exp,env)
  coexistence <- sum(post$feasibility_results)/nrow(post)
  
  post$feasibility_results <- as.factor(post$feasibility_results)
  
  theta_omega <- ggplot(post)+
    geom_point(mapping = aes(x=omega_results,y=theta_results, col= feasibility_results), show.legend = FALSE) +  
    scale_color_manual(values = c(col2,col1)) +
    xlim(0,0.8)  +
    ylim(0,100)   +
    labs(col = "Coexistence") +
    xlab("Omega (niche differences)")+
    ylab("Theta (fitness differences)")+
    annotate(geom="text",label=coexistence, size=6,   x = .25, y = 95)+
    theme_bw()


  
  return(theta_omega)
  
}




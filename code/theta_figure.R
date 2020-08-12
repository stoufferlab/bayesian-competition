#test for figure regarding theta vs. omega
require(gridExtra)
source("code/model_toolbox.R")

make_figure_theta <- function(vero_model,vero_function, vero_exp,trcy_model,trcy_function,trcy_exp,env,title){
 
   col1<-rethinking::col.alpha("mediumseagreen",.3)
   col2<-rethinking::col.alpha("grey50",.3)
  
  post        <-  feasibility_wrapper(vero_model,vero_function, vero_exp,trcy_model,trcy_function,trcy_exp,env)
  coexistence <- sum(post$feasibility_results)/nrow(post) 
  coexistence <- round(coexistence, digits = 2)
  
  
  mean_post <- apply(post,2,median)
  
  post$feasibility_results <- as.factor(post$feasibility_results)
  
  theta_omega <- ggplot(post)+
    geom_point(mapping = aes(x=omega_results,y=theta_results, col= feasibility_results), show.legend = FALSE) +  
    scale_color_manual(values = c(col2,col1)) +
    xlim(0,0.5)  +
    ylim(0,85)   +
    labs(col = "Coexistence") +
    xlab("")+
    ylab("")+
    annotate(geom="text",label=coexistence, size=5,   x = .25, y = 80)+
    ggtitle(title)
    theme_bw()
    

  
  pp<-theta_omega + geom_point(mapping= aes(x= mean_post[1], y=mean_post[3]), col= "goldenrod", size=3)
  
  return(pp)
  
}

#test for figure regarding theta vs. omega
require(gridExtra)
source("code/model_toolbox.R")

make_figure_theta <- function(vero_model,trcy_model,vero_growth,trcy_growth,v1,v2){
  col1<-rethinking::col.alpha("dodgerblue3",1)
  col2<-rethinking::col.alpha("firebrick3",.3)
  
  post        <- posterior_feasibility(vero_model,trcy_model,vero_growth,trcy_growth,v1,v2)
  coexistence <- sum(post$feasibility_results)/nrow(post)
  
  post$feasibility_results <- as.factor(post$feasibility_results)
  fixed_coexistence        <- as.data.frame(fixed_feasibility(vero_model,trcy_model,vero_growth,trcy_growth,v1,v2))
  
  theta_omega <-ggplot(post)+
    geom_point(mapping = aes(x=omega_results,y=theta_results, col= feasibility_results), show.legend = FALSE) +  
    scale_color_manual(values = c(col2,col1)) +
    xlim(0,0.8)  +
    ylim(0,100)   +
    labs(col = "Coexistence") +
    xlab("Omega (niche differences)")+
    ylab("Theta (fitness differences)")+
    annotate(geom="text",label=coexistence, size=6,   x = .25, y = 95)+
    theme_bw()


  theta_omega_final <- theta_omega + geom_point(fixed_coexistence, mapping = aes(x=fixed_coexistence[1,],y=fixed_coexistence[4,]),col="gold", size=2)
  
  return(theta_omega_final)
  
}


BV_BV<-make_figure_theta(BEV_vero,BEV_trcy,vero_growth,trcy_growth,"BV","BV")

RC_RC<-make_figure_theta(RC_vero,RC_trcy,vero_growth,trcy_growth,"RC","RC")

RC_BV<-make_figure_theta(RC_vero,BEV_trcy,vero_growth,trcy_growth,"RC","BV")

BV_RC<-make_figure_theta(BEV_vero,RC_trcy,vero_growth,trcy_growth,"BV","RC")

LA_BV<-make_figure_theta(LAW_vero,BEV_trcy,vero_growth,trcy_growth,"LA","BV")

LA_RC<-make_figure_theta(LAW_vero,RC_trcy,vero_growth,trcy_growth,"LA","RC")

LA_LA<-make_figure_theta(LAW_vero,LAW_trcy,vero_growth,trcy_growth,"LA","LA")

BV_LA<-make_figure_theta(BEV_vero,LAW_trcy,vero_growth,trcy_growth,"BV","LA")

RC_LA<-make_figure_theta(RC_vero,LAW_trcy,vero_growth,trcy_growth,"RC","LA")



pdf(file="omega_more", width = 14)
grid.arrange(BV_BV, RC_BV, LA_BV,
             BV_RC, RC_RC, LA_RC,
             BV_LA, RC_LA, LA_LA, 
             ncol=3, nrow=3)
dev.off()
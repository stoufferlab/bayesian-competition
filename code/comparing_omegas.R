#playing with omegas
library(brms)
library(ggplot2)
library(ggpubr)

#We source everything known to human kind...

source("code/gg_theme.R")
source("code/read_models.R")
source("code/model_toolbox.R")
source("code/model_combo.R")
source("code/feasibility_toolbox.R")


#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033



Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol =  1.88129e-25 )
  #calculate the prob density
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  #out <- log10(d[1]) + n * log10(2)
  # return( d[1]^(1 / n))
  return(d)
}




# the list of models over to iterate
vero_models <- list( vero_bh_multispecies_poisson.rds)

trcy_models<- list(trcy_lv_multispecies_poisson.rds)

model_grid_sunny_linear<- combined_models(vero_models = vero_models,
                                   trcy_models = trcy_models,
                                   si =si,
                                   gi =gi,
                                   gj =gj,
                                   sj=sj,
                                   env=FALSE)




Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

  
  model_grid_sunny_saavedra<- combined_models(vero_models = vero_models,
                                            trcy_models = trcy_models,
                                            si =si,
                                            gi =gi,
                                            gj =gj,
                                            sj=sj,
                                            env=FALSE)
  
  
  Omega <- function(alpha) {
    S <- nrow(alpha)
    Sigma <- solve(t(alpha) %*% alpha)
    # m <- matrix(0, S, 1)
    # a <- matrix(0, S, 1)
    # b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    
    d[1]^(1/S)
  }


  
  
  model_grid_sunny_song<- combined_models(vero_models = vero_models,
                                              trcy_models = trcy_models,
                                              si =si,
                                              gi =gi,
                                              gj =gj,
                                              sj=sj,
                                              env=FALSE)
  
  
col1 <- rethinking::col.alpha("mediumseagreen", .5)
col2 <- rethinking::col.alpha("grey50", .5)

feasibility_plot_sunny<-ggplot(model_grid_sunny_linear) +
  geom_point(
    mapping = aes(
      x = omega_results,
      y = theta_results,
      col = as.factor(feasibility_results)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_mean,
    y = Theta_mean
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  geom_abline(slope = 180, intercept = 0, col="grey50")+
  geom_point(model_grid_sunny_saavedra, mapping = aes(x= omega_results,
                                                      y= theta_results,
                                                      col= as.factor(feasibility_results)), show.legend = FALSE)+
  geom_point(model_grid_sunny_song, mapping = aes(x= omega_results,
                                                      y= theta_results,
                                                      col= as.factor(feasibility_results)), show.legend = FALSE)+
  xlim(-2,.5)+
  ylim(0, 90)+
  facet_grid(trcy_model~vero_model)

feasibility_plot_sunny

#####

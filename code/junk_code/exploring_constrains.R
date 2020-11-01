source("code/feasibility_toolbox.R")
source("code/model_toolbox.R")




#survival and germination for Vero (i) and Trcy(j)
gi<-.372
si<-.556
gj<-.258
sj<-.033



posterior_feasibility<-function(vero_model,trcy_model,si,gi,sj,gj,env){
  
  #for the mean omega and theta we get the  alpha matrix for mean parameter values
  #for an environmental condition
  mean_alpha_matrix <- get_fixed_alphas(vero_model = vero_model,
                                        trcy_model = trcy_model,
                                        gi = gi,
                                        gj = gj,
                                        env = env)
  #as well as r1
  vero_growth <- get_fixed_growth(model = vero_model,
                                  s =si,
                                  g =gi, 
                                  env = env)
  
  #and r2
  trcy_growth <- get_fixed_growth( model = trcy_model,
                                   s = sj,
                                   g = gj,
                                   env =env)
  
  #and now we can calculate the mean feasibility without constraints
  
  Omega_mean <- Omega(mean_alpha_matrix)
  Theta_mean <- theta(mean_alpha_matrix, c(vero_growth, trcy_growth))
  
  #if we constrain it
  constraints<- list(vero_model$constraints, trcy_model$constraints)
  proportion_mean<- calculate_constrained_proportion(alpha = mean_alpha_matrix, constraints = constraints)
  Omega_prime_mean<- Omega_mean * proportion_mean
  # for the posterior feasibility
  trcy_post<-posterior_parameters(model = trcy_model, s = sj,g = gj)
  vero_post<-posterior_parameters(model = vero_model, s = si, g = gi)
  
  num_posterior<- identical(nrow(vero_post),nrow(trcy_post))
  if(num_posterior){
    
    omega_results       <-c()
    feasibility_results <-c()
    growth_results      <-c()
    theta_results       <-c()
    proportion_results <-c()
    omega_prime_results <-c()
    
    for( i in 1:nrow(vero_post)){
      
      #we get the corresponding posterior values, vero first, trcy second, gi (vero), gj(trcy)
      alpha  <- alpha_matrix(vero_row=vero_post[i,],trcy_row =trcy_post[i,],gi=gi,gj=gj,env=env)
      
      if(env){
        r1 <- vero_post$env_growth[i]
        r2 <- trcy_post$env_growth[i]
      }else{
        r1 <- vero_post$growth[i]
        r2 <- trcy_post$growth[i]
      }
      
      
      #And estimate the feasability domain unconstrained
      omega       <-Omega(alpha)
      feasibility <-test_feasibility(alpha,c(r1,r2))
      theta       <-theta(alpha,c(r1,r2))
      
      constraints<- list(vero_model$constraints, trcy_model$constraints)
      #and constrained
      proportion <- calculate_constrained_proportion(alpha = alpha, constraints = constraints)
      omega_prime <- omega * proportion
      # #we save it 
      omega_results       <-c(omega_results,omega)
      feasibility_results <-c(feasibility_results, feasibility)
      theta_results       <-c(theta_results,theta)
      proportion_results  <-c(proportion_results, proportion)
      omega_prime_results <- c(omega_prime_results, omega_prime)
    }
    
    
    pp<- as.data.frame(cbind(omega_results,feasibility_results,theta_results, proportion_results, omega_prime_results))
    pp$Omega_mean <- Omega_mean
    pp$Theta_mean <- Theta_mean
    pp$proportion_mean <- proportion_mean
    pp$Omega_prime_mean <- Omega_prime_mean
    
    pp$vero_model <- vero_model$name
    pp$trcy_model <- trcy_model$name
    return(pp)
  }else{warning("Posterior distributions are not the same length")}
  
  
}

#sensu chuliang
post <- posterior_feasibility(vero_model = vero_bh_multispecies_poisson.rds,
                              trcy_model = trcy_bh_multispecies_poisson.rds,
                              si = si,
                              gi = gi,
                              sj = sj, 
                              gj = gj,
                              env = FALSE)

#sensu alba
post_a <- posterior_feasibility(vero_model = vero_bh_multispecies_poisson.rds,
                              trcy_model = trcy_bh_multispecies_poisson.rds,
                              si = si,
                              gi = gi,
                              sj = sj, 
                              gj = gj,
                              env = FALSE)

plot(post$omega_results, post$omega_prime_results, cex=.5, col= rethinking::col.alpha("blue", 0.2), xlab="Omega", ylab="Constrained Omega")


plot(post_a$omega_results, post_a$omega_prime_results, cex=.5, col= rethinking::col.alpha("blue", 0.2), xlab="Omega", ylab="Constrained Omega")

col1 <- rethinking::col.alpha("mediumseagreen", .05)
col2 <- rethinking::col.alpha("grey50", .05)



ggplot(post_a)+
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
  xlim(0,0.3)+
  ylim(0, 60)

ggplot(post_a)+
  geom_point(
    mapping = aes(
      x = omega_prime_results,
      y = theta_results,
      col = as.factor(feasibility_results)
    ),
    show.legend = FALSE
  ) +
  geom_point( mapping = aes(
    x = Omega_prime_mean,
    y = Theta_mean
  ), col= "goldenrod3") +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  geom_abline(slope = 180, intercept = 0, col="grey50")+
  xlim(0,0.3)+
  ylim(0, 60)





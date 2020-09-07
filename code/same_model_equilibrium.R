library(brms)
library(ggplot)


posterior_pairs_pois <- posterior_parameters(model = vero_bh_pairs_poisson, growth_fun = bev_growth,equilibrium_fun = bh_equilibrium, s = si, g = gi,exp_param = FALSE)

posterior_pairs_zipois <- posterior_parameters(model = vero_bh_pairs_zipoisson, growth_fun = bev_growth,equilibrium_fun = bh_equilibrium, s = si, g = gi,exp_param = FALSE)

posterior_multi_zipois <- posterior_parameters(model = vero_bh_multispecies_zippoisson, growth_fun = bev_growth,equilibrium_fun = bh_equilibrium, s = si, g = gi,exp_param = FALSE)

posterior_multi_pois <- posterior_parameters(model = vero_bh_multispecies_poisson, growth_fun = bev_growth,equilibrium_fun = bh_equilibrium, s = si, g = gi,exp_param = FALSE)


p <- ggplot(posterior_pairs_pois) +
  geom_density(mapping = aes(x=equilibrium), fill="darkblue", alpha=.8 )+
  theme_bw()+ 
  geom_density(mapping = aes(x=equilibrium), fill="firebrick", alpha=.7, data = posterior_pairs_zipois ) +
  geom_density(mapping = aes(x=equilibrium), fill="darkorange", alpha=.7, data = posterior_multi_zipois )+
 geom_density(mapping = aes(x=equilibrium), fill="cyan4", alpha=.8, data = posterior_multi_pois )+
  xlim(0,600)
source("code/model_toolbox.R")
source("code/read_models.R")
source("code/feasibility_toolbox.R")
source("code/gg_theme.R")
#many posteriors
gi<-.372
si<-.556
gj<-.258
sj<-.033
col1 <- rethinking::col.alpha("mediumseagreen", .3)
col2 <- rethinking::col.alpha("grey50", .3)


#for one model of trcy calculate all of the possible model combinations of vero

fixed_row_model <- function(vero_models, trcy_model, si=si,gi = gi,sj = sj,gj = gj,env = FALSE){
  one_row<-lapply(vero_models,function(m){
    post <- feasibility_wrapper(vero_model = m, 
                                trcy_model = trcy_model,
                                si = si,
                                gi = gi,
                                sj = sj,
                                gj = gj,
                                env = env)
    return(post)
  })
  
  
  all_posteriors <- do.call(rbind, one_row)
  return(all_posteriors)

}




combined_models<-function(vero_models, trcy_models, si=si,gi = gi,sj = sj,gj = gj,env = FALSE){
  many_rows <- lapply(trcy_models, function(m2){
    one_row<-fixed_row_model(vero_models = vero_models,
                             trcy_model = m2,
                             si = si,
                             gi = gi,
                             sj = sj,
                             gj = gj,
                             env = env)
    return(one_row)
  })
  
  all_rows<- do.call(rbind, many_rows)
  return(all_rows)
}




test<- combined_models(vero_models = list(vero_bh_multispecies_poisson.rds, 
                                          vero_lv_multispecies_poisson.rds),
                        trcy_models = list(trcy_bh_multispecies_poisson.rds,
                                           trcy_lv_multispecies_poisson.rds),
                       si =si,
                       gi =gi,
                       gj =gj,
                       sj=sj,
                       env=FALSE)







pp<-fixed_row_model(vero_models=list(vero_bh_multispecies_poisson.rds,
                    vero_lv_multispecies_poisson.rds,
                    vero_rc_multispecies_poisson.rds,
                    vero_hs_multispecies_poisson.rds),
                    si = si,
                    trcy_model = trcy_bh_multispecies_poisson.rds, 
                    gi = gi,
                    sj = sj,
                    gj = gj,
                    env = FALSE)



ggplot(pp) +
  geom_point(
    mapping = aes(
      x = omega_results,
      y = theta_results,
      col = as.factor(feasibility_results)
    ),
    show.legend = FALSE
  ) +
  theme_alba +
  scale_color_manual(values = c(col2, col1)) +
  facet_wrap(~vero_model)





ggplot(post_2) +
  geom_point(
    mapping = aes(
      x = omega_results,
      y = theta_results,
      col = as.factor(feasibility_results)
    ),
    show.legend = FALSE
  ) +
  theme_alba +
  scale_color_manual(values = c(col2, col1))

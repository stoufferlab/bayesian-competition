#scripts to itearte over model combinations, to get all posteriors of omega and theta

#for one model of trcy calculate all of the possible model combinations of vero
# vero models is a list of models
fixed_row_model <- function(vero_models,
                            trcy_model,
                            si,
                            gi,
                            Ni_max,
                            sj,
                            gj,
                            Nj_max,
                            env = FALSE,
                            bounded = TRUE){
  one_row<-lapply(vero_models,function(m,
                                       trcy_model,
                                       si,
                                       gi,
                                       sj,
                                       gj,
                                       Ni_max,
                                       Nj_max,
                                       env,
                                       bounded){
    post <- posterior_feasibility(vero_model=m,
                                  trcy_model = trcy_model,
                                  si = si,
                                  gi = gi,
                                  sj = sj,
                                  gj = gj,
                                  Ni_max  = Ni_max,
                                  Nj_max = Nj_max,
                                  env = env,
                                  bounded = bounded)
    
    return(post)
  },  trcy_model = trcy_model,
  si = si,
  gi = gi,
  sj = sj,
  gj = gj,
  Ni_max  = Ni_max,
  Nj_max = Nj_max,
  env = env,
  bounded = bounded )
  
  
  all_posteriors <- do.call(rbind, one_row)
  return(all_posteriors)

}

#trcy models is a list of models as well


combined_models<-function(vero_models,
                          trcy_models,
                          si,
                          gi,
                          Ni_max,
                          sj,
                          gj,
                          Nj_max,
                          env = FALSE,
                          bounded = TRUE){
    many_rows <- lapply(trcy_models, function(m2,
                                              vero_models,
                                              si,
                                              gi,
                                              Ni_max,
                                              sj,
                                              gj,
                                              Nj_max,
                                              env,
                                              bounded){
    one_row<-fixed_row_model(vero_models = vero_models,
                             trcy_model = m2,
                             si=si,
                             gi = gi,
                             Ni_max = Ni_max,
                             sj = sj,
                             gj = gj,
                             Nj_max = Nj_max,
                             env = env,
                             bounded = bounded)
    return(one_row)
  }, vero_models = vero_models,
  si=si,
  gi = gi,
  Ni_max = Ni_max,
  sj = sj,
  gj = gj,
  Nj_max = Nj_max,
  env = env,
  bounded = bounded)
  
  all_rows<- do.call(rbind, many_rows)
  return(all_rows)
}






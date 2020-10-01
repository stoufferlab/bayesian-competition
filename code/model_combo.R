#scripts to itearte over model combinations, to get all posteriors of omega and theta

#for one model of trcy calculate all of the possible model combinations of vero
# vero models is a list of models
fixed_row_model <- function(vero_models, trcy_model, si=si,gi = gi,sj = sj,gj = gj,env = FALSE){
  one_row<-lapply(vero_models,function(m){
    post <- posterior_feasibility(vero_model = m, 
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

#trcy models is a list of models as well


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






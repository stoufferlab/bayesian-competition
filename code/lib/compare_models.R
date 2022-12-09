
require(brms)

#Function to compare models, you give it as many brms objects as you want!
model_compare <- function(..., sort.by = "loo_weights", woody=NULL){
  models <- list(...)
  
  # print a warning if only one model is supplied
  if (length(models) == 1) {
    warning("Only one model is supplied. Comparing one model is not meaningful.")
  }
  
  if(is.null(woody)){
    obs_ind <- NULL
  }else{
    # this assumes all models were fit with the same data!
    data <- models[[1]]$data

    #under which environmental conditions?
    obs_ind <- which(data$woody==woody)
  }

  # switching to lapply from for-loop
  ic <- lapply(models, function(m) {
    w <- loo_subsample(m, observations=obs_ind, loo_approximation='waic')
    wa <- w$estimates["looic", "Estimate"]  # in newer brms access with "estimate", switch back to $waic if you're still using older brms
    
    l <- loo_subsample(m, observations=obs_ind)
    lo <- l$estimates["looic", "Estimate"]
    
    return(c(waic = wa, looic = lo))
  })
  ic <- do.call(rbind, ic)
  
  waic_weights <- exp(-0.5*ic[,"waic"])/sum(exp(-0.5*ic[,"waic"]))
  # model_weights(..., weights = "waic", observations=obs_ind)
  loo_weights <- waic_weights <- exp(-0.5*ic[,"looic"])/sum(exp(-0.5*ic[,"looic"]))
  # model_weights(..., weights = "loo", observations=obs_ind)
  
  out <-
    cbind(waic = ic[, "waic"],
          waic_weights,
          looic = ic[, "looic"],
          loo_weights)

  rownames(out) <- lapply(models, function(x) x$name)
  # out <- out[order(out[, sort.by], decreasing=TRUE),]  # sort.by can be alternatively "looic" now
  return(as.data.frame(out))
  
}

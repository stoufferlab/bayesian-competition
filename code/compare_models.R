#Function to compare models, you give it as many brms objects as you want!




model_compare <- function(..., sort.by = "waic"){
  models <- list(...)
  
  # print a warning if only one model is supplied
  if (length(models) == 1) {
    warning("Only one model is supplied. Comparing one model is not meaningful.")
  }
  
  # switching to lapply from for-loop
  ic <- lapply(models, function(m) {
    w <- waic(m)
    wa <- w$estimates["waic", "Estimate"]  # in newer brms access with "estimate", switch back to $waic if you're still using older brms
    
    l <- loo(m)
    lo <- l$estimates["looic", "Estimate"]
    
    return(c(waic = wa, looic = lo))
  })
  ic <- do.call(rbind, ic)
  
  waic_weights <- model_weights(..., weights = "waic")
  loo_weights <- model_weights(..., weights = "loo")
  
  out <-
    cbind(waic = ic[, "waic"],
          waic_weights,
          looic = ic[, "looic"],
          loo_weights)
  out <- out[order(out[, sort.by]),]  # sort.by can be alternatively "looic" now
  return(out)
  
}
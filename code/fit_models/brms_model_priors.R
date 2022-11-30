
# set prior distributions for brms models

require(brms)

# which priors should be used
model_prior <- function(model_name){
  # for intrinsic fecundity (in log space)
  prior_lambda= c(
    prior(cauchy(0, 5), lb=0, nlpar = "lambdaopen"),
    prior(cauchy(0, 5), lb=0, nlpar = "lambdawoody")
  )

  # for focal-focal interaction ceofficients to be positively constrained
  prior_alphas_pos= c(
    prior(normal(0, 1), lb=0, nlpar = "alphaiiopen"),
    prior(normal(0, 1), lb=0, nlpar = "alphaiiwoody"),
    prior(normal(0, 1), lb=0, nlpar = "alphaijopen"),
    prior(normal(0, 1), lb=0, nlpar = "alphaijwoody"),
    prior(normal(0, 1), nlpar = "alphaikopen"),
    prior(normal(0, 1), nlpar = "alphaikwoody")
  )

  # for focal-focal interaction coefficients to be positive or negative
  prior_alphas= c(
    prior(normal(0, 1), nlpar = "alphaiiopen"),
    prior(normal(0, 1), nlpar = "alphaiiwoody"),
    prior(normal(0, 1), nlpar = "alphaijopen"),
    prior(normal(0, 1), nlpar = "alphaijwoody"),
    prior(normal(0, 1), nlpar = "alphaikopen"),
    prior(normal(0, 1), nlpar = "alphaikwoody")
  )

  if(model_name == 'interaction_free')
    return(prior_lambda)
  else if(model_name == "beverton_holt")
    return(c(prior_lambda, prior_alphas_pos))
  else if(model_name == "ricker")
    return(c(prior_lambda, prior_alphas_pos))
  else
    stop("invalid model_name passed to model_priors()")
}

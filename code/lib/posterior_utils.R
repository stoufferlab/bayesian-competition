
require(brms)
require(posterior)

#Functions to extract parameters from bayesian fits and convert them to equivalent phenomenological parameters

#function to generate a posterior growth rate (and an environmental growth rate) for each point of the posterior, for one model. s and g are the seed survival and germination rates. Model is the brms model fit. It spits out the posterior parameters, the posterior growth and environmental growth and equilibrium abundances and environmental equilibrium abundances
posterior_parameters <- function(
  model,
  s,
  g,
  model_name
){

    # extract samples from the model posterior for all parameters
    post        <- posterior::as_draws_df(model)

    # the lambda parameter is positively constrained during sampling
    lambda_open   <- post$b_lambdaopen_Intercept
    lambda_woody  <- post$b_lambdawoody_Intercept

    # only works when the model has interaction coefficients
    if("alphaiiopen" %in% names(model$formula$pforms)) {
      # intraspecific alphas in control and treatment conditions
      alphaii_open  <- post$b_alphaiiopen_Intercept
      alphaii_woody <- post$b_alphaiiwoody_Intercept

      # interspecific alphas in control and treatment conditions
      alphaij_open  <- post$b_alphaijopen_Intercept
      alphaij_woody <- post$b_alphaijwoody_Intercept
    }else{
      # intraspecific alphas in control and treatment conditions
      alphaii_open  <- NA
      alphaii_woody <- NA

      # interspecific alphas in control and treatment conditions
      alphaij_open  <- NA
      alphaij_woody <- NA
    }
    
    # create a posterior sample of "more interpretable" model parameters
    posterior <-
      as.data.frame(
        cbind(
          lambda_open,
          lambda_woody,
          alphaii_open,
          alphaii_woody,
          alphaij_open,
          alphaij_woody
        )
      )

    # add in the other alpha parameter if it exists (never used but inferred so kept for completeness)
    if("alphaikopen" %in% names(model$formula$pforms)) {
      posterior$alphaik_open  <- post$b_alphaikopen_Intercept
      posterior$alphaik_woody <- post$b_alphaikwoody_Intercept
    }else{
      posterior$alphaik_open  <- NA
      posterior$alphaik_woody <- NA
    }

    # parameter transformations are necessary for the Beverton-Holt model
    if(model_name=="Beverton-Holt"){
      alpha_params <- grep("alpha",colnames(posterior),value=TRUE)
      open_alphas  <- grep("open", alpha_params, value=TRUE)
      woody_alphas <- grep("woody", alpha_params, value=TRUE)
      posterior[,open_alphas] <- sweep(
        posterior[,open_alphas],
        1,
        posterior$lambda_open,
        "/"
      )
      posterior[,woody_alphas] <- sweep(
        posterior[,woody_alphas],
        1,
        posterior$lambda_woody,
        "/"
      )
      posterior$lambda_open  <- 1/posterior$lambda_open
      posterior$lambda_woody <- 1/posterior$lambda_woody
    }

    #These definitions are true for all models
    a_open  <- posterior$lambda_open*g
    a_woody <- posterior$lambda_woody*g
    b       <- 1-((1-g)*s)
    
    # composite eta parameter (seeds produced per seed lost to death or germination)
    posterior$eta_open <- a_open / b
    posterior$eta_woody <- a_woody / b

    #to calculate different growth rates based on the model we define growth based on the name of the model
    if(model_name=="Beverton-Holt"){
      r_open  <- -1 + (a_open/b)
      r_woody <- -1 + (a_woody/b)
    }
    
    if(model_name=="Lotka-Volterra"){
      r_open  <- 1 - (b/a_open) 
      r_woody <- 1 - (b/a_woody)
    }
    
    if(model_name=="Ricker"){
      r_open  <- log(a_open/b)
      r_woody <- log(a_woody/b)
    }

    if(model_name=="Null"){
      r_open <- r_woody <- NA
    }

    # add in the composite "growth rates"
    posterior$r_open  <- r_open
    posterior$r_woody <- r_woody
    
    # add in the germination and seed survival parameters for use in feasibility calculations
    posterior$g <- g
    posterior$s <- s

    return(posterior)
}

# function to return an alpha matrix based on one row of the posterior for each species
# woody is a binary that tells it if to take into consideration the environmental differences
alpha_matrix <- function(
  vero_row,
  trcy_row,
  woody
){
  #if env, then alphas are alpha + alpha_env 
  if(woody){
    alpha11 <- vero_row$alphaii_woody
    alpha21 <- trcy_row$alphaii_woody
    alpha12 <- vero_row$alphaij_woody
    alpha22 <- trcy_row$alphaii_woody
  }else{
    alpha11 <- vero_row$alphaii_open
    alpha21 <- trcy_row$alphaij_open
    alpha12 <- vero_row$alphaij_open
    alpha22 <- trcy_row$alphaii_open
  }
  alpha <- matrix(c(alpha11,alpha21,alpha12,alpha22),ncol=2,nrow=2)
  
  # TODO: should we get rid of g here to keep things consistent with the nomenclature in the paper?
  # alpha   <- sweep(alpha,MARGIN=2,STAT=c(vero_row$g,trcy_row$g),FUN="*")
  return(alpha)  
}

# function to return an alpha matrix based on one row of the posterior for each species
# woody is a binary that tells it if to take into consideration the environmental differences
growth_rates <- function(
  vero_row,
  trcy_row,
  woody
){
  if(woody){
    r1 <- vero_row["r_woody"]
    r2 <- trcy_row["r_woody"]
  }else{
    r1 <- vero_row["r_open"]
    r2 <- trcy_row["r_open"]
  }
  return(as.numeric(c(r1,r2)))
}


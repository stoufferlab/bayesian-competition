#Prepare data set of VERO vs TRCY
#set data and priors for fitting brms models with the environment as the fixed effect. 
library(brms)
source("code/clean_data.R")

#Set the priors all models are going to use 
prior = c(
  prior(normal(0, 10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 10), nlpar = "lambda")
)


prior_exp = c(
  prior(normal(0, 10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 10), nlpar = "lambda"),
  prior(normal(0, 10), nlpar = "beta")
)


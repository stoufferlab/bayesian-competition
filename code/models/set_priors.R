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
# 
# prior_vero_alt = c(
#   prior(normal(0,10), nlpar = "alphaii"),
#   prior(normal(0, 10), nlpar = "alphaij"),
#   prior(normal( 0, 10), nlpar = "lambdai")
# )
# 
# prior_vero_exp = c(
#   prior(normal(0,1), nlpar = "alphaii"),
#   prior(normal(0, 1), nlpar = "alphaij"),
#   prior(normal(0, 10), nlpar = "lambdai"),
#   prior(normal(0, 1), nlpar = "bi")
# )
# 
# prior_trcy = c(
#   prior(normal(0,10), nlpar = "alphajj"),
#   prior(normal(0, 10), nlpar = "alphaji"),
#   prior(normal(0, 10), nlpar = "lambdaj")
# )
# 
# 
# prior_trcy_alt = c(
#   prior(normal(0,10), nlpar = "alphajj"),
#   prior(normal(0, 10), nlpar = "alphaji"),
#   prior(normal(trcy_mean, 50), nlpar = "lambdaj")
# )
# 
# 
# prior_trcy_exp = c(
#   prior(normal(0,10), nlpar = "alphajj"),
#   prior(normal(0, 10), nlpar = "alphaji"),
#   prior(normal(0, 10), nlpar = "lambdaj"),
#   prior(normal(0, 10), nlpar = "bj")
# )
# 
# 

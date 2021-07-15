#Prepare data set of VERO vs TRCY
#set data and priors for fitting brms models with the environment as the fixed effect. 
library(brms)

# #Set the priors all models are going to use 
# prior_pairs = c(
#   prior(normal(0, 10), nlpar = "alphaii"),
#   prior(normal(0, 10), nlpar = "alphaij"),
#   prior(normal(0, 10), nlpar = "lambda")
# )

prior_multi= c(
  prior(normal(0, 0.1), nlpar = "alphaii"),
  prior(normal(0, 0.1), nlpar = "alphaij"),
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, 0.1), nlpar = "alphaik")
)


prior_lotka= c(
  prior(normal(0, 0.01), nlpar = "alphaii"),
  prior(normal(0, 0.01), nlpar = "alphaij"),
  prior(normal(0, 0.1), nlpar = "lambda"),
  prior(normal(0, 0.01), nlpar = "alphaik")
)


# 
# prior_exp_pairs = c(
#   prior(normal(0, 10), nlpar = "alphaii"),
#   prior(normal(0, 10), nlpar = "alphaij"),
#   prior(normal(0, 10), nlpar = "lambda"),
#   prior(normal(0, .1 ), nlpar = "beta")
# )
# 


 prior_exp_multi = c(
   prior(normal(0, 10), nlpar = "alphaii"),
   prior(normal(0, 10), nlpar = "alphaij"),
   prior(normal(0, 10), nlpar = "lambda"),
  prior(normal(0,  0.1), nlpar = "beta"),
   prior(normal(0, 10), nlpar = "alphaik")
 )


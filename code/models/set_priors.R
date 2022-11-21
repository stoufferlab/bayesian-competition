
#set prior distributions for brms models
library(brms)

prior_multi= c(
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, 0.1), nlpar = "alphaii"),
  prior(normal(0, 0.1), nlpar = "alphaij"),  
  prior(normal(0, 0.1), nlpar = "alphaik")
)

# prior_exp_multi = c(
#   prior(normal(0, 10), nlpar = "lambda"),
#   prior(normal(0, 10), nlpar = "alphaii"),
#   prior(normal(0, 10), nlpar = "alphaij"),
#   prior(normal(0, 10), nlpar = "alphaik")
#   # prior(normal(0,  0.1), nlpar = "beta")
# )

# prior_lotka= c(
#   prior(normal(0, 0.01), nlpar = "alphaii"),
#   prior(normal(0, 0.01), nlpar = "alphaij"),
#   prior(normal(0, 0.1), nlpar = "lambda"),
#   prior(normal(0, 0.01), nlpar = "alphaik")
# )

# 
# prior_exp_pairs = c(
#   prior(normal(0, 10), nlpar = "alphaii"),
#   prior(normal(0, 10), nlpar = "alphaij"),
#   prior(normal(0, 10), nlpar = "lambda"),
#   prior(normal(0, .1 ), nlpar = "beta")
# )
# 

#Prepare data set of VERO vs TRCY
#set data and priors for fitting brms models with the environment as the fixed effect. 
library(brms)
source("code/clean_data.R")
#This script works with vero and trcy, so we remove plantago
rm(vero_plde)
#vero = i
#trcy = j
vero_focal<- vero_trcy %>% filter(focal == "V")
trcy_focal<- vero_trcy %>% filter(focal == "T")

#remove the original data set so I do not get confused later
rm(vero_trcy)

#Set the priors all models are going to use 
prior_vero = c(
  prior(normal(0,10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 10), nlpar = "lambdai")
)

prior_vero_exp = c(
  prior(normal(0,10), nlpar = "alphaii"),
  prior(normal(0, 10), nlpar = "alphaij"),
  prior(normal(0, 10), nlpar = "lambdai"),
  prior(normal(0, 10), nlpar = "bi")
)




prior_trcy = c(
  prior(normal(0,100), nlpar = "alphajj"),
  prior(normal(0, 100), nlpar = "alphaji"),
  prior(normal(0, 100), nlpar = "lambdaj")
)




prior_trcy_exp = c(
  prior(normal(0,10), nlpar = "alphajj"),
  prior(normal(0, 10), nlpar = "alphaji"),
  prior(normal(0, 10), nlpar = "lambdaj"),
  prior(normal(0, 10), nlpar = "bj")
)



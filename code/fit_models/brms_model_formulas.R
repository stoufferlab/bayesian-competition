
# brms nonlinear model structures for phenomenological models of density-dependent fecundity

require(brms)

# what non-linear model formula should be used
model_formula <- function(model_name){
  # INTERACTION-FREE MODEL
  interaction_free<-bf(
    totalseeds ~ (1-woody)*lambdaopen + woody*lambdawoody,
    lambdaopen ~ 1,
    lambdawoody ~ 1,
    nl=TRUE
  )

  # BEVERTON-HOLT MODEL
  beverton_holt<-bf(
    totalseeds ~ (((1-woody)*lambdaopen + woody*lambdawoody) + (((1-woody)*alphaiiopen + woody*alphaiiwoody)*conspecifics) + (((1-woody)*alphaijopen + woody*alphaijwoody)*heterospecifics) + (((1-woody)*alphaikopen + woody*alphaikwoody)*totalother))^(-1),
    lambdaopen  ~ 1,
    lambdawoody  ~ 1,
    alphaiiopen ~ 1,
    alphaijopen ~ 1,
    alphaikopen ~ 1,
    alphaiiwoody ~ 1,
    alphaijwoody ~ 1,
    alphaikwoody ~ 1,
    nl=TRUE
  )

  # RICKER MODEL
  ricker<- bf(
    totalseeds ~ ((1-woody)*lambdaopen + woody*lambdawoody) * exp( -(((1-woody)*alphaiiopen + woody*alphaiiwoody)*conspecifics) - (((1-woody)*alphaijopen + woody*alphaijwoody)*heterospecifics) - (((1-woody)*alphaikopen + woody*alphaikwoody)*totalother)),
    lambdaopen  ~ 1,
    lambdawoody  ~ 1,
    alphaiiopen ~ 1,
    alphaijopen ~ 1,
    alphaikopen ~ 1,
    alphaiiwoody ~ 1,
    alphaijwoody ~ 1,
    alphaikwoody ~ 1,
    nl=TRUE
  )

  if(model_name == 'interaction_free')
    return(interaction_free)
  else if(model_name == "beverton_holt")
    return(beverton_holt)
  else if(model_name == "ricker")
    return(ricker)
  else
    stop("invalid model_name passed to model_formula()")
}
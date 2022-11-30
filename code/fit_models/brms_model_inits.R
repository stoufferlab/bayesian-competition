
# get a good guess for initial parameter values via a glm fit

# require(manyglm)

# which priors should be used
model_init <- function(model_name, data, chains=1){
  glm_fit <- glm(
    totalseeds ~ 1 + woody,
    data=data,
    family=poisson
  )
  if(model_name == "interaction_free"){
    init <- list(
      b_lambdaopen = array(data=exp(coef(glm_fit)["(Intercept)"])),
      b_lambdawoody = array(data=exp(coef(glm_fit)["(Intercept)"] + coef(glm_fit)["woody"]))
    )
    return(lapply(1:chains, function(x) init))
  }else if(model_name == "ricker"){
    mformula <- as.formula(totalseeds ~ 1 + woody + conspecifics * woody + heterospecifics * woody + totalother * woody)
    lfunc <- 'log'
    start <- c(
      (coef(glm_fit)["(Intercept)"]),
      (coef(glm_fit)["(Intercept)"] + coef(glm_fit)["woody"]),
      rep(0,6)
    )
  }else if(model_name == "beverton_holt"){
    mformula <- as.formula(totalseeds ~ 1 + woody + conspecifics * woody + heterospecifics * woody + totalother * woody)
    lfunc <- 'inverse'
    start <- c(
      1/exp(coef(glm_fit)["(Intercept)"]),
      1/exp(coef(glm_fit)["(Intercept)"] + coef(glm_fit)["woody"]),
      rep(0,6)
    )
  }else{
    stop("invalid model_name passed to model_inits()")    
  }

  # fit a glm to the data
  glm_fit <- glm(
    mformula,
    data=data,
    family=poisson(link=lfunc),
    start=start
  )

  # return appropriately specified initial values
  if(model_name == "beverton_holt"){
    init <- list(
      b_lambdaopen = array(
        data=coef(glm_fit)["(Intercept)"]
      ),
      b_lambdawoody = array(
        data=(coef(glm_fit)["(Intercept)"] + coef(glm_fit)["woody"])
      ),
      b_alphaiiopen = array(
        data=0.001 #coef(glm_fit)["conspecifics"]
      ),
      b_alphaiiwoody = array(
        data=0.001 #(coef(glm_fit)["conspecifics"]+coef(glm_fit)["woody:conspecifics"])
      ),
      b_alphaijopen = array(
        data=0.001 #coef(glm_fit)["heterospecifics"]
      ),
      b_alphaiiwoody = array(
        data=0.001 #(coef(glm_fit)["heterospecifics"]+coef(glm_fit)["woody:heterospecifics"])
      ),
      b_alphaikopen = array(
        data=0.001 #coef(glm_fit)["totalother"]
      ),
      b_alphaikwoody = array(
        data=0.001 #(coef(glm_fit)["totalother"]+coef(glm_fit)["woody:totalother"])
      )
    )
  }else if(model_name == "ricker"){
    init <- list(
      b_lambdaopen = array(
        data=exp(coef(glm_fit)["(Intercept)"])
      ),
      b_lambdawoody = array(
        data=exp(coef(glm_fit)["(Intercept)"] + coef(glm_fit)["woody"])
      ),
      b_alphaiiopen = array(
        data=0.001 #*coef(glm_fit)["conspecifics"]
      ),
      b_alphaiiwoody = array(
        data=0.001 #*(coef(glm_fit)["conspecifics"]+coef(glm_fit)["woody:conspecifics"])
      ),
      b_alphaijopen = array(
        data=0.001 #*coef(glm_fit)["heterospecifics"]
      ),
      b_alphaiiwoody = array(
        data=0.001 #*(coef(glm_fit)["heterospecifics"]+coef(glm_fit)["woody:heterospecifics"])
      ),
      b_alphaikopen = array(
        data=0.001 #*coef(glm_fit)["totalother"]
      ),
      b_alphaikwoody = array(
        data=0.001 #*(coef(glm_fit)["totalother"]+coef(glm_fit)["woody:totalother"])
      )
    )
  }

  return(lapply(1:chains, function(x) init))
}
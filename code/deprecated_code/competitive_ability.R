###competitive ability of different competition models, following the definition of Hart et al. 

bev_competitive_ability<-function(s,g,lambda, alphaii, alphaij ){
  a <- lambda*g
  b <- 1 - ((1-g)*s)
  growth <- (a/b) - 1
  competition <- sqrt(alphaii *alphaij  )
  ability     <- growth/competition
  return(competition)
}

lotka_competitive_ability<-function(s,g,lambda, alphaii, alphaij){
  a<- lambda*g
  b<- 1-( (1-g)*s )
  growth <- 1  -(b/a) 
  competition <- sqrt(alphaii *alphaij  )
  ability <- growth/competition
  return(competition)
}

ricker_competitive_ability<- function(s,g,lambda,alphaii,alphaij){
  a<- lambda*g
  b<- 1- ( (1-g)*s)
  growth<- log(a/b)
  competition <- sqrt(alphaii *alphaij  )
  ability <- growth/competition
  return(competition)
}


law_competitive_ability<- function(s,g,lambda,b_param,alphaii,alphaij){
  a<- lambda*g
  b<- 1- ( (1-g)*s)
  growth<- -1 + ( (a/b)^ (1/b_param))
  competition <- sqrt(alphaii *alphaij  )
  ability <- growth/competition
  return(competition)
}

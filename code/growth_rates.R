#Code to figure out the limits of a growth rate
# the term lambda is not multiplied by g because its estimated value already includes this
bev_growth<-function(s,g,lambda){
  a <- lambda*g
  b <- 1 - ((1-g)*s)
  r <- (a/b) - 1
  return(r)
}
    
lotka_growth<-function(s,g,lambda){
  a<- lambda*g
  b<- 1-( (1-g)*s )
  r<- 1  -(b/a) 
  return(r)
}


ricker_growth<-function(s,g,lambda){
  a<- lambda*g
  b<- 1- ( (1-g)*s)
  r<- log(a/b)
  return(r) 
}

#because this model takes 2 posterior distributions..
law_growth<-function(s,g,lambda,b_param){
  a<- lambda*g
  b<- 1- ( (1-g)*s)
  r<- -1 + (a/b)^b_param
  return(r) 
}


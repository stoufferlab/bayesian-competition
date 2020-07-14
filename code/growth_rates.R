#Code to figure out the limits of a growth rate

bev_growth<-function(s,g,lambda){
  a <- g*lambda
  b <- 1 - ((1-g)*s)
  r <- (a/b) - 1
  return(r)
}
    
lotka_growth<-function(s,g,lambda){
  b<- 1-( (1-g)*s  )
  a<- g*lambda
  r<- -(b/a) + 1
  return(r)
}


ricker_growth<-function(s,g,lambda){
  b<- 1- ( (1-g)*s)
  a<- g*lambda
  r<- -log(b/a)
  return(r) 
}



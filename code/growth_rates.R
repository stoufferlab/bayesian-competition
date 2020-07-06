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
#s goes from 0 -1
#g goes from 0 - 1
#lambda goes from 0, inf
bev_growth(1,0,10)
lotka_growth(.5,1,1000)
ricker_growth(1,1,.1)




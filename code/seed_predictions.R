gi<-.372
si<-.556
gj<-.258
sj<-.033

Ni <- seq(0,100,1)
Nj <- seq(0,100,1)



BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
BEV_trcy <- readRDS("~/bayesian-competition/BEV_trcy.RDS")

RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
RC_trcy <- readRDS("~/bayesian-competition/RC_trcy.RDS")

LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")
LAW_trcy <- readRDS("~/bayesian-competition/LAW_trcy.RDS")



fixed_model<-function(model){
  model_coef<-fixef(model)
  coef<-as.matrix(model_coef[,1])
  coef<-t(coef)
  params<-as.data.frame(coef)
  return(params)
}


vero_bev <- fixed_model(BEV_vero)
vero_law <- fixed_model(LAW_vero)
#
# bev_seeds<-function(si,gi,gj,lambda,alphaii,alphaij,Ni,Nj){
#   seed_bank   <- (( 1-gi)*si)
#   fecundity   <- lambda
#   competition <- 1 - (alphaii*Ni) - (alphaij*Nj)
#   prod        <- (fecundity/competition)
#   seeds       <- seed_bank + (fecundity/competition)
#   return (prod)
# }



bev_seeds<-function(lambda,alphaii,alphaij,Ni,Nj){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  fecundity   <- lambda *gi
  competition <- 1 + (alphaii*gi*Ni) + (alphaij*gj*Nj)
  seeds       <-  (fecundity/competition)

  return (seeds)
}




law_seeds<-function(lambda,alphaii,alphaij,Ni,Nj,b){
  gi<-.372
  si<-.556
  gj<-.258
  sj<-.033
  
  fecundity   <- lambda *gi
  competition <- (1 + (alphaii*gi*Ni) + (alphaij*gj*Nj) ) ^b
  seeds       <-  (fecundity/competition)
  
  return (seeds)
}








v1 <- bev_seeds( lambda=vero_bev$lambda_Intercept,alphaii =  vero_bev$alphaii_Intercept,
                alphaij = vero_bev$alphaij_Intercept,Ni = Ni, Nj = 0)
v2 <- bev_seeds( lambda=vero_bev$lambda_Intercept,alphaii =  vero_bev$alphaii_Intercept,
                 alphaij = vero_bev$alphaij_Intercept,Ni = 0, Nj = Nj)




v3 <- bev_seeds( lambda=vero_law$lambda_Intercept,alphaii =  vero_law$alphaii_Intercept,
                 alphaij = vero_law$alphaij_Intercept,Ni = Ni, Nj = 0)
v4 <- bev_seeds( lambda=vero_law$lambda_Intercept,alphaii =  vero_law$alphaii_Intercept,
                 alphaij = vero_law$alphaij_Intercept,Ni = 0, Nj = Nj)






plot(Ni,v1, ylim=c(0,5),type="l", lwd=2, xlab="Neighbors", ylab="seeds")
lines(Nj,v2, lwd=2, lty=2)
lines(Ni,v3, lwd=2, lty=1, col= "dodgerblue4")







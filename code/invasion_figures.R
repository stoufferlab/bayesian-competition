
ff<-coexistence_condition(3)





colour_generator<-function(red,green,blue ,alpha){
  r<-red/255
  g<-green/255
  b<-blue/255
  col<-rgb(r,g,b,alpha)
}

col1<-colour_generator(52,118,208,.5)
col2<-colour_generator(208,142,52,.5)



coex_point<-function(lambdai,alphaii,alphaij,lambdaj,alphajj,alphaji,col){
  fitness_ratio<-(lambdai-1)/(lambdaj-1)  * sqrt(  (alphaji*alphajj) / (alphaij*alphaii) )
  niche_overlap<-sqrt(  (alphaij*alphaji)  / (alphaii*alphaij))
  points(niche_overlap,fitness_ratio,col=col,pch=20)
  ff<-c(niche_overlap, fitness_ratio)
}



niche_overlap<-function(lambdai,alphaii,alphaij,lambdaj,alphajj,alphaji,col){
  fitness_ratio<-(lambdai-1)/(lambdaj-1)  * sqrt(  (alphaji*alphajj) / (alphaij*alphaii) )
  niche_overlap<-sqrt(( alphaji/alphaii) * (alphaij/alphajj)  )
  #niche_overlap<-sqrt(  (alphaij*alphaji)  / (alphaii*alphaij))
  #points(niche_overlap,fitness_ratio,col=col,pch=20)
  #ff<-c(niche_overlap, fitness_ratio)
  points(niche_overlap,lambdai,col=col,pch=20)
}




plot(NULL,xlim=c(0,2),ylim=c(0,1),xlab="Niche overlap", ylab="Growth rate of sp.i", main="")

parnames1<-c("b_lambda_Intercept","b_alphaii_Intercept","b_alphaij_Intercept")
parnames2<-c("b_lambda_Intercept","b_alphajj_Intercept", "b_alphaji_Intercept")
posterior_points_alphas<-function(model1,parnames1,model2,parnames2, n_samples){
  sub<-seq(1,n_samples,1)
  post1<-posterior_samples(model1,parnames1,exact_match = T, subset=sub)
  lambdai<-post1[,1]
  alphaii<-post1[,2]
  alphaij<-post1[,3]
  post2<-posterior_samples(model2,parnames2,exact_match = T, subset=sub)
  lambdaj<-post2[,1]
  alphajj<-post2[,2]
  alphaji<-post2[,3]
  points_refill<-data.frame("lambdai"=lambdai, "alphaii"=alphaii,"alphaij"=alphaij,"lambdaj"=lambdaj , "alphajj"=alphajj, "alphaji"=alphaji )
  points_refill
}

 posterior<-posterior_points_alphas(bev_vero,parnames1,bev_trcy,parnames2,5000)
 posterior2<-posterior_points_alphas(exp_vero,parnames1,exp_trcy,parnames2,5000)
#  mcmc_scatter(posterior,pars=c("alphaii","alphajj"))
#  posterior_2<-posterior_points_alphas(exp_vero,parnames1,exp_trcy,parnames2,100)
#  
#  coex_point()
# add_coex_points<-function(posterior,col){
#   nn<-nrow(posterior)
#   for(i in 1:nn){
#     coex_point(posterior[,1],posterior[,2], posterior[,3],posterior[,4],posterior[,5],posterior[,6],col=col)
#   }
# }



add_niche_points<-function(posterior,col){
  nn<-nrow(posterior)
  for(i in 1:nn){
    niche_overlap(posterior[i,1],posterior[i,2], posterior[i,3],posterior[i,4],posterior[i,5],posterior[i,6],col=col)
  }
}




 add_niche_points(posterior,col1)
 add_niche_points(posterior2, col2)
 
 
 print(summary(exp_vero),digits=4)
 print(summary(exp_trcy),digits=4)
 
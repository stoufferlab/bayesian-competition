require(tidyverse)
#Consider that we estimate the parameters of both species with a B-H model
#we have the following interaction matrix
#Sp.1 facilitates species two but limits itself, sp.2 limits itself and species 1.
alpha<-matrix(c(1,-0.7, 0.7, 1), nrow = 2, ncol = 2)
#We can calculate the feasibility domain by sampling growth rates, asuming no "constraints".
#Lets also assume that Infinity is a big number like 1e6

no_constraints<- list(c(-1e6,1e6),c(-1e6,1e6))
Nsample <- 1e5
# We sample N number of growth rates
r <- no_constraints %>% 
  map(~runif(Nsample, min = .[1], max = .[2])) %>% 
  unlist() %>% 
  matrix(ncol = Nsample, byrow = TRUE)


#And test and save which ones are feasible given our alpha
test_feasibility <- function(alpha,r){
  out <- all(solve(alpha,r)>0)
  return(out)
}

feasibility <- c()
r1 <- c()
r2 <- c()
for(i in 1:Nsample){

  feasibility[i] <-test_feasibility(alpha,r[,i])
  r1[i]<-r[,i][1]
  r2[i]<-r[,i][2]
  
}

#Originally you would calculate Omega by:
omega<- mean(feasibility)

#But we know that is an over estimation of the feasibility domain, because there are some constraints on the possible values of the growth rates.
#So we can sample inside the feasibility domain to see which proportion of it is compatible with our constraints. First we define the feasibility domaain

full_domain<-cbind(feasibility,r1,r2) %>% as.data.frame()
feasibility_domain <- filter(full_domain, feasibility ==1)

#Sp.1 and 2 can only have growth rates that go from -1, to Inf. So lets take out any combination that is outside those limits
feasibility_domain_constrained<- filter(feasibility_domain, r1 >= -1 & r2 >= -1 )

#And we can estimate now Omega prime as the new proportion of those values that fall inside the feasibility domain AND are within our constraints
omega_prime<-sum(feasibility_domain_constrained$feasibility) / Nsample


#This way we guarantee that Omega prime is the same size or smaller than Omega.
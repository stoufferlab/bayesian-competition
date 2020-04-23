

colour_generator<-function(red,green,blue ,alpha){
    r<-red/255
    g<-green/255
    b<-blue/255
    col<-rgb(r,g,b,alpha)
}

col1<-colour_generator(0,0,0,.3)



coexistence<-function(niche,max_fitness){
    fitness_ratio<-seq(0,max_fitness,.01)
    inverse_niche<-1/niche
    coexist<-c()
    for(i in 1:length(fitness_ratio)){
        
        if(fitness_ratio[i]>niche & fitness_ratio[i]<inverse_niche){
            coexist[i]<-fitness_ratio[i]
        }else{
            coexist[i]<-0
        }
    }
    coexist
}


coexistence_condition<-function(max_fitness){
    niche_overlap<-seq(0,1,.01)
    ff1<-c()
    ff2<-c()
    n1<-c()
    n2<-c()
    aa<-c()
    bb<-c()
    
    plot(NULL, xlim=c(0,1), ylim=c(0,3), xlab="Niche Overlap",ylab="Fitness ratio")
    for(i in 1:length(niche_overlap)){
        coex<-coexistence(niche_overlap[i],max_fitness)
        fitness_conditions<-coex[which(!coex==0)]
        fmin<-min(fitness_conditions)
        ff1<-c(ff1,fmin)
        n1<-c(n1,niche_overlap[i])
        fmax<-max(fitness_conditions)
        ff2<-c(ff2,fmax)
        fits<-c(fmin,fmax)
        niches<-c(niche_overlap[i], niche_overlap[i])
        aa<-c(aa,niches)
        bb<-c(bb,fits)
       
     
    }
  # lines(n1,ff1,lwd=1.5)
  # lines(n1,ff2,lwd=1.5)
   f1<-ff1[1:length(ff1)-1]
   f2<-ff2[1:length(ff2)-1]
   ni<-n1[1:length(n1)-1]
   polygon(c(ni,rev(ni)), c(f2,rev(f1)),col="gray79" )
}


fig<-coexistence_condition(3)



    
#the plde and vero data set
dat1<-read.csv("data/PLDEvVero.FinalDatasetFromAIP.csv")
#the trycy data set
dat2<-read.csv("data/VEROvTRCY.FinalDatasetFromAIP.csv")




vero_plde<-data.frame(focal=dat1$Focal,label=dat1$label, verotreatment=dat1$Vtrt, veroplanted=dat1$VEROseeds,verodensity=dat1$Vn, pldetreatment=dat1$Ptrt, pldeplanted=dat1$PLDEseeds, pldedensity=dat1$Pn, totalother=dat1$total.other, originalenv=dat1$Original.Env, exposedenv=dat1$Exp.Env, replicate=dat1$replicate, totalseeds= dat1$Seeds.total)


vero_trcy<-data.frame(focal=dat2$focal,label=dat2$label, verotreatment=dat2$Nv, veroplanted=dat2$Veroseeds.planted, verodensity=dat2$Nv.1,  trcytreatment=dat2$Nt, trcyplanted=dat2$TRCYseeds.planted, trcydensity=dat2$Nt.1, totalother=dat2$totalother, originalenv=dat2$Original.Env, exposedenv=dat2$Exp.Env, replicate=dat2$replicate, totalseeds=dat2$number.of.seeds)

rm(dat1,dat2)

# used for the color scheme
library(here)
# library(viridis)

testing <- FALSE

feasibility_posteriors <- read.csv(here('data/results/combined_feasibility_posteriors.csv'))

# # remove undetected feasibility domains
# feasibility_posteriors <- subset(feasibility_posteriors, domain_detected)

feasibility_posteriors$woody <- ifelse(
    feasibility_posteriors$woody,
    "Open",
    "Woody"
)

# # outcome table
# outcome_table <- table(
#     feasibility_posteriors$biological_outcome,
#     paste(feasibility_posteriors$vero_model, feasibility_posteriors$trcy_model, feasibility_posteriors$woody, sep="_")
# )

# # order the table
# outcome_table <- outcome_table[c("i_too_big","i","BOTH",'j','j_too_big'),]

# # normalize each column
# outcome_probs <- sweep(
#     outcome_table,
#     2,
#     colSums(outcome_table),
#     "/"
# )

# lets make a lovely figure of all this too while we're at it
if(!testing){
    pdf(
        file=here('figures/posterior_structural_stability.pdf'),
        width=14,
        height=10
    )
}else{

}

layout.matrix <- matrix(
    c(1,2,3,4),
    nrow = 2,
    ncol = 2,
    byrow = T
)
layout(mat = layout.matrix, heights=c(1,1), widths = c(1,1))
# layout(mat = c(1), heights=c(2), widths = c(2))

par(oma = c(2, 4, 4, 5))

par(mar = c(2, 5, 2, 5))

for(tmodel in c("Beverton-Holt","Ricker")){
    for(vmodel in c("Beverton-Holt","Ricker")){
        if(tmodel=="Ricker"){
            ann<-TRUE
        }else{
            ann<-FALSE
        }
        dd.tmp <- subset(feasibility_posteriors, vero_model==vmodel & trcy_model==tmodel)
        dd.tmp$woody <- as.factor(dd.tmp$woody)
        library(vioplot)
        vioplot(
            log(area_feasible) ~ woody,
            dd.tmp,
            cex.axis=2.1,
            ann=ann,
            ylab=''
            # ylim=c(-2,14)
            # xlab=c('Open','Woody')
            # ann=TRUE
        )
        
        # mylevels <- levels(dd.tmp$woody)
        # for(i in 1:length(mylevels)){
        #     thislevel <- mylevels[i]
        #     thisvalues <- log(dd.tmp$area_feasible[dd.tmp$woody==thislevel])
        #     myjitter <- jitter(rep(i, length(thisvalues)), amount=0.15)
        #     points(myjitter, thisvalues)
        # }
   
        # small_outcome_probs <- outcome_probs[,c(paste0(vmodel,"_",tmodel,"_1"),paste0(vmodel,"_",tmodel,"_0"))]
        # barplot(
        #     small_outcome_probs,
        #     horiz=TRUE,
        #     col=viridis(5),
        #     names.arg = c("Woody","Open"),
        #     cex.names=1.75,
        #     cex.axis=1.5
        # )


        # # model labels
        if(vmodel=="Beverton-Holt"){
            # mtext(tmodel,side=2,cex=2,line=5.5)
            # title(ylab='Size of the biologically constrained feasibility domain',cex.lab=1.5,xpd=TRUE)
            if(tmodel=="Ricker"){
                text(5.65,3,tmodel,adj=c(0.5,0.5),cex=2.2,xpd=NA,srt=270)
                text(
                    x=0,
                    y=14,
                    label=expression('Size of the biologically constrained feasibility domain, log('*italic(A[beta])*')'),
                    # side=2,
                    cex=2.5,
                    # line=3.5,
                    srt=90,
                    xpd=NA
                )
            }else{
                text(5.65,4,tmodel,adj=c(0.5,0.5),cex=2.2,xpd=NA,srt=270)
            }
        }
        if(tmodel=="Beverton-Holt"){
            if(vmodel=="Beverton-Holt"){
                text(1.5,15,vmodel,adj=c(0.5,0.5),cex=2.2,xpd=NA)
            }else{
                text(1.5,10.9,vmodel,adj=c(0.5,0.5),cex=2.2,xpd=NA)
            }
            # mtext(vmodel,side=3,cex=2,line=1)
        }else{

        }

        # # add a symbol to indicate median prediction
        # median_outcome_0 <- subset(feasibility_posteriors, median & vero_model == vmodel & trcy_model == tmodel & woody==0)$biological_outcome
        # median_outcome_1 <- subset(feasibility_posteriors, median & vero_model == vmodel & trcy_model == tmodel & woody==1)$biological_outcome
        # # stop()
        # midpoints_0 <- (2*cumsum(small_outcome_probs[,2])-small_outcome_probs[,2])/2
        # midpoints_1 <- (2*cumsum(small_outcome_probs[,1])-small_outcome_probs[,1])/2
        # points(midpoints_0[median_outcome_0],1.9,bg='white',pch=21,cex=5)
        # points(midpoints_1[median_outcome_1],0.7,bg='white',pch=21,cex=5)

        # stop()
    }
}
# legend(
#     "topright",
#     fill=viridis(5),
#     c(
#         expression(italic("G. rosea")*"*"),
#         expression(italic("G. rosea")),
#         "Coexistence",
#         expression(italic("T. cyanopetala")),
#         expression(italic("T. cyanopetala")*"*")
#     ),
#     pt.bg=viridis(5),
#     xpd=NA,
#     inset=c(-0.55,-0.22),
#     # bty='n',
#     cex=1.5,
#     title="Predicted outcome"
# )
text(
    -0,
    27.5,
    expression(italic("Goodenia rosea")),
    xpd=NA,
    # side=3,
    # outer=TRUE,
    cex=2.9,
    # line=1.5
)
par(xpd=TRUE)
text(
    3,
    9,
    expression(italic("Trachymene cyanopetala")),
    xpd=NA,
    adj=c(0.5,0.5),
    cex=2.9,
    srt=270
)

if(!testing){
    dev.off()
}

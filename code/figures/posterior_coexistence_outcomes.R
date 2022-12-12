
# used for the color scheme
library(here)
library(viridis)

feasibility_posteriors <- read.csv(here('data/results/combined_feasibility_posteriors.csv'))

# outcome table
outcome_table <- table(
    feasibility_posteriors$biological_outcome,
    paste(feasibility_posteriors$vero_model, feasibility_posteriors$trcy_model, feasibility_posteriors$woody, sep="_")
)

# order the table
outcome_table <- outcome_table[c("i_too_big","i","BOTH",'j','j_too_big'),]

# normalize each column
outcome_table <- sweep(
    outcome_table,
    2,
    colSums(outcome_table),
    "/"
)

# lets make a lovely figure of all this too while we're at it
pdf(
    file=here('figures/posterior_coexistence_outcomes.pdf'),
    width=14,
    height=10
)

layout.matrix <- matrix(
    c(1,2,3,4),
    nrow = 2,
    ncol = 2,
    byrow = T
)
layout(mat = layout.matrix, heights=c(1,1), widths = c(2,2))

par(oma = c(2, 7, 7, 10))

par(mar = c(2, 3, 2, 5))

for(tmodel in c("Beverton-Holt","Ricker")){
    for(vmodel in c("Beverton-Holt","Ricker")){
        small_outcome_table <- outcome_table[,c(paste0(vmodel,"_",tmodel,"_1"),paste0(vmodel,"_",tmodel,"_0"))]
        barplot(
            small_outcome_table,
            horiz=TRUE,
            col=viridis(5),
            names.arg = c("Woody","Open"),
            cex.names=1.75,
            cex.axis=1.5
        )
        # model labels
        if(vmodel=="Beverton-Holt")
            mtext(tmodel,side=2,cex=2,line=3.5)
        if(tmodel=="Beverton-Holt")
            mtext(vmodel,side=3,cex=2)
        # add a symbol to indicate median prediction
        median_outcome_0 <- subset(feasibility_posteriors, median & vero_model == vmodel & trcy_model == tmodel & woody==0)$biological_outcome
        median_outcome_1 <- subset(feasibility_posteriors, median & vero_model == vmodel & trcy_model == tmodel & woody==1)$biological_outcome
        # stop()
        midpoints_0 <- (2*cumsum(small_outcome_table[,2])-small_outcome_table[,2])/2
        midpoints_1 <- (2*cumsum(small_outcome_table[,1])-small_outcome_table[,1])/2
        points(midpoints_0[median_outcome_0],1.9,bg='white',pch=21,cex=5)
        points(midpoints_1[median_outcome_1],0.7,bg='white',pch=21,cex=5)
    }
}
legend(
    "topright",
    fill=viridis(5),
    c(
        expression(italic("G. rosea")*"*"),
        expression(italic("G. rosea")),
        "Coexistence",
        expression(italic("T. cyanopetala")),
        expression(italic("T. cyanopetala")*"*")
    ),
    pt.bg=viridis(5),
    xpd=NA,
    inset=c(-0.55,-0.22),
    # bty='n',
    cex=1.5,
    title="Predicted outcome"
)
mtext(expression(italic("Trachymene cyanopetala")),side=2,outer=TRUE,cex=2.5,line=4)
mtext(expression(italic("Goodenia rosea")),side=3,outer=TRUE,cex=2.5,line=1.5)

dev.off()

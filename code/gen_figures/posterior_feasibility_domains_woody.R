
# packages used for parallelization
library(doMC)
registerDoMC(cores = 7)

# define functions used below
require(here)
source(here("code/lib/feasibility_utils.R"))
source(here("code/lib/posterior_utils.R"))

# general code that finds and loads all candidate model fits
source(here("code/lib/load_models.R"))

# parameters not estimated via inference
source(here("data/germination/extra_params.R"))

# tuning parameters that mostly control speed of these calculations
set.seed(123454)
subsample <- 50
desired_feasible <- 5000
max_samples <- 50000

# extract derived posteriors of pop dyn model parameter values
# replace those posteriors with the medians
posteriors <- list(vero=list(),trcy=list())
for(model in c('ricker','beverton_holt')){

    vero_model <- focal_models[["vero"]][[model]]
    posteriors[["vero"]][[model]] <- posterior_parameters(
        vero_model,
        s=seed_survival["vero"],
        g=germination["vero"],
        model_name=vero_model$name
    )
    trcy_model <- focal_models[["trcy"]][[model]]
    posteriors[["trcy"]][[model]] <- posterior_parameters(
        trcy_model,
        s=seed_survival["trcy"],
        g=germination["trcy"],
        model_name=trcy_model$name
    )
}

only.median <- FALSE
for(model in c('ricker','beverton_holt')){
    # attach the median to the front of the posteriors
    posteriors[["vero"]][[model]] <- rbind(
        apply(posteriors[["vero"]][[model]], 2, median),
        posteriors[["vero"]][[model]]
    )
    posteriors[["trcy"]][[model]] <- rbind(
        apply(posteriors[["trcy"]][[model]], 2, median),
        posteriors[["trcy"]][[model]]
    )
    if(only.median){
        posteriors[["vero"]][[model]] <- posteriors[["vero"]][[model]][1,,drop=FALSE]
        posteriors[["trcy"]][[model]] <- posteriors[["trcy"]][[model]][1,,drop=FALSE]
    }else{
        sub_samp <- sample.int(nrow(posteriors[["vero"]][[model]])-1, subsample) + 1
        posteriors[["vero"]][[model]] <- posteriors[["vero"]][[model]][c(1,sub_samp),]
        posteriors[["trcy"]][[model]] <- posteriors[["trcy"]][[model]][c(1,sub_samp),]
    }
}

for(woody in c(1)){

filename <- paste0(
    "figures/posterior_feasibility_domains_",
    ifelse(woody,"woody","open"),
    ".pdf"
)
pdf(
    file=here(filename),
    width=7,
    height=7
)
layout.matrix <- matrix(
    c(1,2,3,4),
    nrow = 2,
    ncol = 2,
    byrow = T
)
layout(mat = layout.matrix, heights=c(1,1), widths = c(1,1))

par(oma = c(2, 2, 0, 0))

par(mar = c(2, 3, 2, 1.5))


xlim <- ylim <- c(0,500)

xmin <- list()
ymin <- list()
xminr <- list()
yminr <- list()
xmax <- list()
ymax <- list()
xmaxr <- list()
ymaxr <- list()
cntr <- 0
for(tmodel in c('beverton_holt','ricker')){
    for(vmodel in c('beverton_holt','ricker')){
        
                model_combo <- paste('v',vmodel, 't',tmodel, woody, sep="_")
                # xmin[[vmodel]] <- Inf
                # ymin[[tmodel]] <- Inf
                # xmax[[vmodel]] <- 0
                # ymax[[tmodel]] <- 0

                # manually tweaked limits based on previously collected xmax, ymax data
                min_inv <- Inf
                if(vmodel=="beverton_holt"){
                    xlim <- c(1/min_inv,140)
                }else if(vmodel=="ricker"){
                    xlim <- c(1/min_inv,70)
                }
                if(tmodel=="beverton_holt"){
                    ylim <- c(1/min_inv,125)
                }else if(tmodel=="ricker"){
                    ylim <- c(1/min_inv,50)
                }

                # message(model_combo)
                # plot scaffold
                plot(
                    x=NA,
                    y=NA,
                    xlim=xlim,
                    ylim=ylim,
                    type='n',
                    xlab="",
                    ylab="",
                    # log='xy',
                    # main=model_combo,
                    cex.lab=1.5
                )
                # dashed lines to guide the eyes
                abline(h=0,lty='dashed',lwd=1.5)
                abline(v=0,lty='dashed',lwd=1.5)

                # to store whether or not the posterior of r is inside or outside the hull
                r_in_hull <- numeric(nrow(posteriors[["vero"]][[vmodel]]))
                for(i in rev(seq.int(nrow(posteriors[["vero"]][[vmodel]])))){
                    alpha  <- alpha_matrix(
                        vero_row = posteriors[["vero"]][[vmodel]][i,],
                        trcy_row = posteriors[["trcy"]][[tmodel]][i,],
                        woody = woody
                    )
                    R_max <- maximum_radius(
                        alpha = alpha, 
                        gN_max = gN_max
                    )
                    # growth rate constraints
                    # TODO: convert this to a function?
                    rconstraints <- list(
                        lower = c(vero=focal_models[["vero"]][[vmodel]]$constraints[1], trcy=focal_models[["trcy"]][[tmodel]]$constraints[1]),
                        upper = c(vero=focal_models[["vero"]][[vmodel]]$constraints[2], trcy=focal_models[["trcy"]][[tmodel]]$constraints[2])
                    )
                    results <- integrate_area(
                        alpha = alpha,
                        gN_max = gN_max,
                        rconstraints = rconstraints,
                        npts = max_samples
                    )$coords
                    # model_combo <- paste(vmodel,tmodel,woody,sep='_')
            
                    bcfd_hull <- grDevices::chull(results)
                    bcfd_hull <- c(bcfd_hull, bcfd_hull[1])
                    bcfd_hull <- results[bcfd_hull,]

                    xmin[[vmodel]] <- min(xmin[[vmodel]], min(bcfd_hull[,1]))
                    ymin[[tmodel]] <- min(ymin[[tmodel]], min(bcfd_hull[,2]))
                    xmax[[vmodel]] <- max(xmax[[vmodel]], max(bcfd_hull[,1]))
                    ymax[[tmodel]] <- max(ymax[[tmodel]], max(bcfd_hull[,2]))

                    r <- growth_rates(
                        vero_row = posteriors[["vero"]][[vmodel]][i,],
                        trcy_row = posteriors[["trcy"]][[tmodel]][i,],
                        woody = woody
                    )

                    xminr[[vmodel]] <- min(xminr[[vmodel]],r[1])
                    yminr[[tmodel]] <- min(yminr[[tmodel]],r[2])
                    xmaxr[[vmodel]] <- max(xmaxr[[vmodel]],r[1])
                    ymaxr[[tmodel]] <- max(ymaxr[[tmodel]],r[2])

                    if(i>1){
                        polygon(
                            bcfd_hull,
                            border = NA,
                            lty = 0,
                            col = scales::alpha("mediumseagreen", alpha = 0.1)
                        )
                    }else{
                        polygon(
                            bcfd_hull,
                            border = "black",
                            lty=1
                            # col = scales::alpha("mediumseagreen", alpha = 0.1)
                        )
                    }

                    r_in_hull_test <- sp::point.in.polygon(
                        r[1],
                        r[2],
                        bcfd_hull$ri,
                        bcfd_hull$rj
                    )
                    r_in_hull[i] <- (r_in_hull_test>0)
                    
                }
                # to plot the vital rates r
                for(i in rev(seq.int(nrow(posteriors[["vero"]][[vmodel]])))){
                    r <- growth_rates(
                        vero_row = posteriors[["vero"]][[vmodel]][i,],
                        trcy_row = posteriors[["trcy"]][[tmodel]][i,],
                        woody = woody
                    )
                    pch <- ifelse(r_in_hull[i],21,23)
                    if(i==1){
                        points(
                            r[1],
                            r[2],
                            pch=pch,
                            col='black',
                            bg='mediumseagreen'
                        )
                    }else{
                        pcol <- ifelse(r_in_hull[i],'white','grey')
                        points(
                            r[1],
                            r[2],
                            pch=pch,
                            col='black', #scales::alpha("black", alpha = 0.1),
                            bg=pcol #scales::alpha("white", alpha = 0.1)
                        )
                    }
                }

                vtext <- switch(
                    vmodel,
                    ricker="Ricker",
                    beverton_holt="Beverton-Holt"
                )
                ttext <- switch(
                    tmodel,
                    ricker="Ricker",
                    beverton_holt="Beverton-Holt"
                )

                # label the axes
                if(cntr %in% c(2,3)){
                    mtext(
                        bquote(.(vtext) ~ "vital rate," ~ italic(nu)["G. rosea"]),
                        side = 1,
                        # outer = TRUE,
                        line=2.75,
                        cex=1.2
                    )
                }
                if(cntr %in% c(0,2)){
                    mtext(
                        bquote(.(ttext) ~ "vital rate," ~ italic(nu)["T. cyanopetala"]),
                        side = 2,
                        line=2.25,
                        # outer = TRUE,
                        cex=1.2
                    ) #, adj = 1)
                }

                cntr <- cntr + 1
    }
}

dev.off()

}
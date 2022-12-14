
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
set.seed(22380)
subsample <- NULL
desired_feasible <- 5000
max_samples <- 50000

# get combined feasibility metrics for all models except the null
vero_models_to_mix <- setdiff(names(focal_models[["vero"]]), c("interaction_free"))
trcy_models_to_mix <- setdiff(names(focal_models[["trcy"]]), c("interaction_free"))

# mix models together and calculate the posteriors
feasibility_posteriors <- foreach(vero_model_name = vero_models_to_mix, .combine = rbind) %:%
    foreach(trcy_model_name = trcy_models_to_mix, .combine=rbind) %dopar% {
    # select the models of interest    
    vero_model <- focal_models[["vero"]][[vero_model_name]]
    trcy_model <- focal_models[["trcy"]][[trcy_model_name]]

    # extract derived posteriors of pop dyn model parameter values
    posteriors <- list()
    posteriors[["vero"]] <- posterior_parameters(
        vero_model,
        s=seed_survival["vero"],
        g=germination["vero"],
        model_name=vero_model$name
    )
    posteriors[["trcy"]] <- posterior_parameters(
        trcy_model,
        s=seed_survival["trcy"],
        g=germination["trcy"],
        model_name=trcy_model$name
    )

    # make sure both posteriors have the same number of samples otherwise we cannot pair them one to one
    if(nrow(posteriors[["vero"]]) != nrow(posteriors[["trcy"]]))
        stop('these posteriors cannot be matched')

    # add in a column to quickly identify the median values
    posteriors[["vero"]]$median <- FALSE
    posteriors[["trcy"]]$median <- FALSE

    # attach medians to front of posteriors
    posteriors[["vero"]] <- rbind(
        apply(posteriors[["vero"]], 2, median),
        posteriors[["vero"]]
    )
    posteriors[["vero"]]$median[1] <- TRUE
    posteriors[["trcy"]] <- rbind(
        apply(posteriors[["trcy"]], 2, median),
        posteriors[["trcy"]]
    )
    posteriors[["trcy"]]$median[1] <- TRUE

    # when !is.null(subsample) use a subsample of non-median values
    if(!is.null(subsample)){
        sub_samp <- sample.int(nrow(posteriors[[1]])-1, subsample) + 1
        posteriors[["vero"]] <- posteriors[["vero"]][c(1,sub_samp),]
        posteriors[["trcy"]] <- posteriors[["trcy"]][c(1,sub_samp),]
    }

    # growth rate constraints
    # TODO: convert this to a function?
    rconstraints <- list(
        lower = c(vero=vero_model$constraints[1], trcy=trcy_model$constraints[1]),
        upper = c(vero=vero_model$constraints[2], trcy=trcy_model$constraints[2])
    )

    environmental_posteriors <- data.frame()
    for(woody in c(0,1)){
        # calculate a bunch of metrics about these two species together
        feasibility_stats <- do.call(rbind,lapply(
            seq.int(nrow(posteriors[["vero"]])),
            function(i){
                message(paste(vero_model$name, trcy_model$name, "woody", woody, "sample", i))
                alpha  <- alpha_matrix(
                    vero_row = posteriors[["vero"]][i, ],
                    trcy_row = posteriors[["trcy"]][i, ],
                    woody = woody
                )
                r <- growth_rates(
                    vero_row = posteriors[["vero"]][i, ],
                    trcy_row = posteriors[["trcy"]][i, ],
                    woody = woody
                )
                R_max <- maximum_radius(
                    alpha = alpha, 
                    gN_max = gN_max
                )
                results <- structural_stability_wrapper(
                    r = r,
                    alpha = alpha,
                    R_max = R_max,
                    gN_max = gN_max,
                    rconstraints = rconstraints,
                    desired_feasible = desired_feasible,
                    max_samples = max_samples
                )
                results <- cbind(
                    data.frame(
                        i="vero",
                        j="trcy",
                        alphaii=alpha[1,1],
                        alphaij=alpha[1,2],
                        alphaji=alpha[2,1],
                        alphajj=alpha[2,2],
                        ri=r[1],
                        rj=r[2],
                        R_max=R_max,
                        giNi_max=gN_max[1],
                        gjNj_max=gN_max[2],
                        median=ifelse(i==1,TRUE,FALSE)
                    ),
                    results
                )
                return(results)
            }
        ))
        feasibility_stats$woody <- woody
        feasibility_stats$vero_model <- vero_model$name
        feasibility_stats$trcy_model <- trcy_model$name
        environmental_posteriors <- rbind(environmental_posteriors, feasibility_stats)
    }
    return(environmental_posteriors)
}

# make a table of predicted coexistence outcomes
feasibility_posteriors$biological_outcome <- NA
for(i in seq.int(nrow(feasibility_posteriors))){
    if(feasibility_posteriors$coexist[i]=="i" && feasibility_posteriors$giNi_mono[i] > feasibility_posteriors$giNi_max[i]){
        feasibility_posteriors$biological_outcome[i] <- "i_too_big"
    }else if(feasibility_posteriors$coexist[i]=="j" && feasibility_posteriors$gjNj_mono[i] > feasibility_posteriors$gjNj_max[i]){
        feasibility_posteriors$biological_outcome[i] <- "j_too_big"
    }else{
        feasibility_posteriors$biological_outcome[i] <- feasibility_posteriors$coexist[i]
    }
}

# write this to a file so that we don't have to wait to regenerate it every bloody time!
write.csv(
    feasibility_posteriors,
    file=here('data/results/combined_feasibility_posteriors.csv'),
    quote=FALSE,
    row.names=FALSE
)

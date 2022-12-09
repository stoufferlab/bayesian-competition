
require(here)

# find all saved model objects
files <-list.files(
	path = here('code/fit_models/model_objects/'),
	pattern = ".rds"
)

# read all model objects into a common list
all_models<-sapply(
	files,
	function(x){
		readRDS(here(paste0('code/fit_models/model_objects/',x)))
	},
	simplify=FALSE
)

# use filenames to split by focal species
focal_models <- list()
focal_models[["vero"]] <- all_models[grepl("vero",names(all_models))]
focal_models[["trcy"]] <- all_models[grepl("trcy",names(all_models))]

# strip the cruft from the model names for ease of use later
for(focal in names(focal_models)){
	names(focal_models[[focal]]) <- gsub(
		paste0(focal,"_"),
		'',
		gsub(
			'.rds','',names(focal_models[[focal]])
		),
	)
}

# add some utility variables regarding the fits
for(focal in names(focal_models)){
	for(mname in names(focal_models[[focal]])){
		if(grepl("beverton",mname)){
			focal_models[[focal]][[mname]]$name <- "Beverton-Holt"
			focal_models[[focal]][[mname]]$constraints <- c(-1, Inf)
		}else if(grepl("ricker",mname)){
			focal_models[[focal]][[mname]]$name <- "Ricker"
			focal_models[[focal]][[mname]]$constraints <- c(-Inf, Inf)
		}else if(grepl("lotka",mname)){
			focal_models[[focal]][[mname]]$name <- "Lotka-Volterra"
			focal_models[[focal]][[mname]]$constraints <- c(-Inf, 1)
		}else if(mname=="interaction_free"){
			focal_models[[focal]][[mname]]$name <- "Null"
			focal_models[[focal]][[mname]]$constraints <- c(NA, NA)
		}
	}
}

# # attach model objects to the environment
# list2env(all_models , envir = .GlobalEnv)

# we remove objects that are no longer required
rm(all_models, files, focal, mname)

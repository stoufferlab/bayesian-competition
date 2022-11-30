
require(tidyverse)
require(here)

# read in the raw vero and trcy data set
dat2<-read.csv(here("data/fecundity/VEROvTRCY.FinalDatasetFromAIP.csv"))

# convert to a tibble
dat2<-as_tibble(dat2)

#rename in a way that is easy to understand
vero_trcy <- dat2 %>%
	rename(
		verotreatment = Nv,
		veroplanted = Veroseeds.planted,
		verodensity = Nv.1,
		trcytreatment = Nt,
		trcyplanted = TRCYseeds.planted,
		trcydensity = Nt.1,
		originalenv = Original.Env,
		exposedenv = Exp.Env,
		totalseeds = number.of.seeds
	)

#make the exposed environment a 0/1 dummy variable and drop NA seed counts
vero_trcy <- vero_trcy %>%
	mutate(
		woody = as.integer(exposedenv=="woody")
	) %>%
	drop_na()

#separate by species
vero_focal <- vero_trcy %>%
	filter(focal == "V") %>%
	rename(
		conspecifics = verodensity,
		heterospecifics = trcydensity
	)

trcy_focal <- vero_trcy %>%
	filter(focal == "T") %>%
	rename(
		conspecifics = trcydensity,
		heterospecifics = verodensity
	)

# remove the intermediate data tables
rm(dat2, vero_trcy)

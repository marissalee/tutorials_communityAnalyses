
require(vegan)
require(dplyr) # this has handy dataframe manipulation tools
require(tidyr) # ditto

# -------------------------------------------------------------------#
# load wood endophyte data

# this dataset has been modified and simplified... the original version has not yet been published
exdata <- readRDS(file = "data/woodEndophyte_exData.RData")

# here's the OTU matrix
otus <- exdata[['otus']]
head(otus)
# the OTU matrix has been minimally processed
# it has been quality checked to remove failed samples (i.e. abundance of reads below a threshold) and remove OTUs that clearly are not fungal
dim(otus)
# there are 50 samples (rows) and 1000 OTUS (cols)

# here's the covariate dataframe (with sample metadata too)
covariates <- exdata[['covariates']]
head(covariates)

# Trim out rare OTUs because they don't contribute much information
minPerc <- 20
numSamps <- dim(otus)[1]
minSamps <- floor( numSamps * (minPerc/100) )
x <- apply(otus > 0, 2, sum) # remove OTUs that don't meet this criteria
otus <- otus[, x > minSamps]
dim(otus)
# now we have 55 OTUs

# select environmental variables you'd like to use to explain OTU composition/abundances
covariates %>%
  select(size, waterperc, density) -> envVars

# -------------------------------------------------------------------#

# Use distance-based redundancy analysis (vegan::capscale) and model selection (ordistep) to determine which wood traits best explain variation in OTU composition

# transform the OTU matrix to de-emphasize taxa with low abundances
otus.t <- decostand(otus, 'hellinger')

# build the unconstrained model
mod0 <- capscale(otus.t ~ 1, data = envVars)

# build the fully constrained model
modfull <- capscale(otus.t ~ ., data=envVars, distance='bray')

# model selection
modstep <- ordistep(mod0, scope=formula(modfull))
modstep
# the constrained model explain ~13% of variation in the data (inertia proportion constrained)

# summary table for the best model
df <- data.frame(term=row.names(modstep$anova), modstep$anova)
df %>%
  separate(term, into=c("drop","term")) %>%
  select(-drop) %>%
  rename('pval'=`Pr..F.`) -> df
df
# model selection suggests that waterperc and size are the only traits that we need in the model to explain OTU composition

# Methods used to analyze associations between wood endophyte OTUs and environmental wood traits
# Part 1: Distance and model-based methods

require(vegan)
require(mvabund)
require(dplyr) # this has handy dataframe manipulation tools
require(tidyr) # ditto

# -------------------------------------------------------------------#
# load data

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

# -------------------------------------------------------------------#

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


# -------------------------------------------------------------------#

# Use GLMs for multivariate abundance data (mvabund::manyglm) and model selection () to determine which wood traits best explain variation in OTU abundances

# specify data for mvabund
fung <- list(abund = otus, x = envVars)
Dat <- mvabund(fung$abund, row.names=row.names(fung$abund)) # so R knows to treat coverDat as multivariate abundance
size <- factor(fung$x$size)
waterperc <- fung$x$waterperc
density <- fung$x$density

# specify mvabund models
# default family = "negative binomial", assumes a quadratic mean-variance relationship and a log-linear relationship between the response variables and any continuous variables
ft.base <- manyglm(Dat ~ size)
ft.m.density <- manyglm(Dat ~ size + density)
ft.m.waterperc <- manyglm(Dat ~ size + waterperc)
mod.list <- list(base = ft.base, 
                 density = ft.m.density, 
                 waterperc = ft.m.waterperc)

# mean-variance plot
meanvar.plot(Dat, xlab="Mean OTU abundance", ylab="Variance in OTU abundance") # with log scale

# check residuals
# each color is associated with an OTU (and there aren't enough colors)
plot(ft.base, which=1, cex=0.5, caption = "", xlab = "")
plot(ft.m.density, which=1, cex=0.5, caption = "", xlab = "")
plot(ft.m.waterperc, which=1, cex=0.5, caption = "", xlab = "")

# compare the AIC value of candidate model to that of the base model

#calculate the total AIC for each model
aic.list <- lapply(mod.list, function(x){
  sum(AIC(x)) # sum AIC values for each OTUId
})

#calculate deltaAIC = base AIC - candidate AIC
baseAIC <- aic.list[['base']]
aic.df <- data.frame(modelName = names(unlist(aic.list)), 
                     AIC = unlist(aic.list))
aic.df %>%
  mutate(deltaAIC = baseAIC - AIC) -> aic.df
aic.df

# lower AIC indicates a better model fit
# so, adding waterperc to the model increases the model fit (decreases the AIC) whereas adding density decreases the model fit (increases the AIC)






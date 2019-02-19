
### Use mvabund and model selection to determine which wood traits best explain variation in OTU abundances
#- Author: Marissa Lee
#- Type: modified project code
#- Test data source: https://github.com/Zanne-Lab/woodEndophytes

# -------------------------------------------------------------------#
# A bit of intro into using mvabund...
# this is largely pinched from http://environmentalcomputing.net/introduction-to-mvabund/
# also see David Warton's mvabund video that was paired with his article in Methods in Ecology and Evolution (http://eco-stats.blogspot.com/2012/03/introducing-mvabund-package-and-why.html)

## How does this method differ from other multivariate analyses? 
# Many commonly used analyses for multivariate data sets (e.g. PERMANOVA, ANOSIM, CCA, RDA etc.) are “distance-based analyses”. 
# This means the first step of the analysis is to calculate a measure of similarity between each pair of samples, thus converting a multivariate dataset into a univariate one.

## Why is this a problem?
# 1) Low statisical power, except for variables with high variance. This means that for variables which are less variable, the analyses are less likely to detect a treatment effect. 
# 2) Does not account for a very important property of multivariate data, which is the mean-variance relationship. Typically, in multivariate datasets like species-abundance data sets, counts for rare species will have many zeros with little variance, and the higher counts for more abundant species will be more variable.

## What is the solution?
# Include an assumption of a mean-variance relationship. 
# Fit a single generalised linear model (GLM) to each response variable (or all response variables with site as a row effect if composition = TRUE) with a common set of predictor variables. 
# We can then use resampling to test for significant community level or species level responses to our predictors.
# Plus, model-based framework makes it easier to check our assumptions and interpret uncertainty around our findings

# -------------------------------------------------------------------#

# load libraries
require(mvabund)
require(dplyr) # this has handy dataframe manipulation tools
require(tidyr) # ditto
require(ggplot2)

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
# Check out how the raw data are distributed

df <- data.frame(samp = row.names(otus), otus)
df %>%
  gather(key = "OTU", value = "num.reads", -samp) -> tmp

ggplot(tmp, aes(x = OTU, y = num.reads)) +
  geom_point() +
  coord_flip()
  
# It looks like some OTUs are much more abundant and variable than others
# It’s probably a good idea to check our mean-variance relationship then... first let's put this data into mvabund-readable format

# -------------------------------------------------------------------#
# Put data in mvabund-friendly format
fung <- list(abund = otus, x = envVars)
Dat <- mvabund(fung$abund, row.names=row.names(fung$abund)) # so R knows to treat coverDat as multivariate abundance
size <- factor(fung$x$size)
waterperc <- fung$x$waterperc
density <- fung$x$density

# -------------------------------------------------------------------#
# Check out the mean-variance relationship
meanvar.plot(Dat, xlab="Mean OTU reads", ylab="Variance in OTU reads") # with log scale

# -------------------------------------------------------------------#
# Specify mvabund models

## Data structure considerations:

# Mean and variance are positively correlated in our data...
# We can deal with this relationship by choosing a family of GLMs with an appropriate mean-variance assumption
# default family = "negative binomial", assumes a quadratic mean-variance relationship and a log-linear relationship between the response variables and any continuous variables
# composition = TRUE/ FALSE (default) fits a separate model to each species. TRUE fits a single model to all variables, including site as a row effect, such that all other terms model relative abundance (compositional effects).

## Biological questions:

# (1) Are endophyte compositions shaped by each of the following host habitat characteristics: 
# stem diamter size (categorical)
# host species density (continuous)
# water content (continuous)

# (2) Assuming that size is a key factor, which is more important is shaping the community -- density or water content?


# -------------------------------------------------------------------#
# Approach for Q1: 

# We can test the multivariate hypothesis of whether species composition varied across the habitats by using the anova function. 
# This gives an analysis of deviance table where we use likelihood ratio tests and resampled p values to look for a significant effect of Habitat on the community data.

ft.size <- manyglm(Dat ~ size, family="negative.binomial", composition = T)
ft.density <- manyglm(Dat ~ density, family="negative.binomial", composition = T)
ft.waterperc <- manyglm(Dat ~ waterperc, family="negative.binomial", composition = T)

# check out what is the in the fitted model object
str(ft.size)
ft.size$coefficients
ft.size$stderr.coefficients

# check residuals
# if you have a linear relationship or fan shape, you might be modeling the wrong distribution
# each color is associated with an OTU (and there aren't enough colors)
plot(ft.size, which=1, cex=0.5, caption = "", xlab = "")
plot(ft.density, which=1, cex=0.5, caption = "", xlab = "")
plot(ft.waterperc, which=1, cex=0.5, caption = "", xlab = "")

# compute the analysis of deviance table -- be careful this can take a while
# resamp = pit.trap (default); Method of resamping the distribution that accounts for correlation in testing
# test = LR (default); Likelihood-Ratio test
# p.uni = none (default); Calculate univariate (e.g. each OTU) test stats and their p-values
# nBoot = 999 (default); Number of bootstrap iterations

an.size <- anova(ft.size, resamp = "pit.trap", test = "LR", p.uni="none", 
      nBoot = 99, show.time = "all") 
an.size

#We can see from this table that there is a significant effect of size (LRT = 153.1, P = 0.01), 
# meaning that the OTU composition clearly differs between the small and large wood stems.

#To examine this further, and see which OTUs are more likely to be found on which size stems, 
# we can run univariate tests for each OTU separately...
# This is done by using the p.uni="adjusted" argument in the anova function. 
# The “adjusted” part of the argument refers to the resampling method used to compute the p values, taking into account the correlation between the response variables.

an.size.uni <- anova(ft.size, resamp = "pit.trap", test = "LR", p.uni="adjusted", 
                 nBoot = 99, show.time = "all") 
an.size.uni
# many individual OTUs are not sensitive to stem size, but a few are like ITSall_OTUa_7145
# now let's figure out whether ITSall_OTUa_7145 increases/decreases with stem size by looking the coefficient estimate

ft.size$coefficients[,'ITSall_OTUa_7145']

# plot the real data to confirm
df <- data.frame(samp = row.names(otus), otus)
df %>%
  gather(key = "OTU", value = "num.reads", -samp) %>%
  filter(OTU == "ITSall_OTUa_7145") %>%
  rename('seq_sampName'='samp') %>%
  left_join(covariates) -> tmp
ggplot(tmp, aes(x = size, y = num.reads)) +
  geom_jitter(width = .2)


# -------------------------------------------------------------------#
# Approach for Q2: Compare the AIC value of candidate model to that of the base model

ft.base <- manyglm(Dat ~ size, family="negative.binomial", composition = T)
ft.m.density <- manyglm(Dat ~ size + density, family="negative.binomial", composition = T)
ft.m.waterperc <- manyglm(Dat ~ size + waterperc, family="negative.binomial", composition = T)
mod.list <- list(base = ft.base, 
                 density = ft.m.density, 
                 waterperc = ft.m.waterperc)

# check residuals
# each color is associated with an OTU (and there aren't enough colors)
plot(ft.base, which=1, cex=0.5, caption = "", xlab = "")
plot(ft.m.density, which=1, cex=0.5, caption = "", xlab = "")
plot(ft.m.waterperc, which=1, cex=0.5, caption = "", xlab = "")

#calculate the total AIC for each model
aic.list <- lapply(mod.list, function(x){
  sum(AIC(x)) # sum AIC values for each OTUId
})
aic.list

#calculate deltaAIC = base AIC - candidate AIC
baseAIC <- aic.list[['base']]
aic.df <- data.frame(modelName = names(unlist(aic.list)), 
                     AIC = unlist(aic.list))
aic.df %>%
  mutate(deltaAIC_base_candidate = baseAIC - AIC) -> aic.df
aic.df

# lower AIC indicates a better model fit
# so, adding waterperc to the model increases the model fit (decreases the AIC) whereas adding density decreases the model fit (increases the AIC)






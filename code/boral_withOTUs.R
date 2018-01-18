
# boral example
# particularly designed for studies with very diverse communities, e.g. OTU data
# custom functions to display ...
# ... mcmc trace plots
# ... distribution of shared enviroment and residual correlations
# ... strength and direction of shared enviroment and residual correlations

require(boral)
require(dplyr) # this has handy dataframe manipulation tools
require(tidyr) # ditto
require(ggplot2) #plotting
require(circlize) #chord diagrams
source("code/boral_withOTUs_auxilaryfunctions.R")

# -------------------------------------------------------------------#
# load data --- this is the same as in part 1, just copied over

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

# Trim out rare OTUs because they don't contribute much information --- this is the same as in part 1, just copied over
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

# Fit model-based constrained and unconstrained ordinations using the boral package

# recode any categorical predictors as a dummy variable
#   in this case, it's size which just has 2 levels
dummy <- data.frame(model.matrix( ~ size - 1, data = envVars))
envVars %>%
  mutate(sizesmall = dummy$sizesmall) %>%
  select(-size) -> envVars

# specify row effects in a rowEffs dataframe
# these are typically blocking variables; ie variables that you have good reason to believe will influence OTU abundances and want to account for it
#   in this case, I'm including 2 row effects...
#   one for each sample (1:50) 
#   --- this is really important for OTU data because we know there is sample-to-sample variation in OTU abundances during lab processing that we need to account for
#   and one for each site from which the samples were collected (1:3)
rowEffs<-data.frame(
  sample = as.numeric(seq(1:nrow(otus))), 
  site = covariates$site
)
#factors need to be coded as numeric
uniqSITE<-unique(rowEffs$site)
for(i in 1:length(uniqSITE)){
  rowEffs$site<-as.character(rowEffs$site)
  rowEffs[rowEffs$site==uniqSITE[i],"site"]<-i
}
rowEffs$site <- as.numeric(rowEffs$site)
row.ids <- matrix(c(rowEffs$sample, rowEffs$site), ncol = 2)

# set mcmc controls
mcmc.controls <- list(n.burnin = 50,
                      n.iteration = 100,
                      n.thin = 1,
                      seed = 1)
# these values are set so that everything runs relatively quickly
# real values are more like...
#n.burnin = 10000
#n.iteration = 100000
#n.thin = 100
# also, its good practice to do model fitting multiple times (9-12?) with new seed values to make sure you are not stuck in a weird part of parameter space

# set priors
prior.controls = list(type = c("normal","normal","normal","uniform"), 
                      hypparams = c(100, 20, 100, 20))

# Fit the model with only latent variables included
# this is analogous to a unconstrained ordination
fit.lv <- boral(y = otus,
             family = "negative.binomial",
             num.lv = 2,
             row.eff = "fixed",
             row.ids = row.ids,
             calc.ics = TRUE, # for model selection
             save.model = TRUE, # to access mcmc samples
             mcmc.control = mcmc.controls,
             prior.control = prior.controls)

# Fit the model with X variation and latent variables included
# this is analogous to a constrained ordination
fit.envlv <- boral(y = otus,
                   X = envVars,
                   family = "negative.binomial",
                   num.lv = 2,
                   row.eff = "fixed",
                   row.ids = row.ids,
                   calc.ics = TRUE, # for model selection
                   save.model = TRUE, # to access mcmc samples
                   mcmc.control = mcmc.controls,
                   prior.control = prior.controls)


# -------------------------------------------------------------------#

# Figure out if your models converged
# note that since these models were run using quick mcmc controls, I am pretty durn sure that they didn't converge
# I'm just going to look at fit.envlv in this section for simplicity but you'll also want to check the convergence of fit.lv and any other models too

# Plot a trace of the MCMC chain for each parameter

# extract the mcmc object (remember, this only gets saved if save.model = TRUE)
mcmc.obj <- as.mcmc.list(fit.envlv$jags.model$BUGSoutput)

# remember that there are typically a ton of parameters estimated in a single boral model
dim(mcmc.obj[[1]]) # this model includes 55 OTUs, 3 covariates, 2 latent variables, and 2 row effects which brings us to 539 parameters

# typically, you can use some pre-packaged functions from within coda to plot the traces
require(coda)
#plot(mcmc.obj) # this will take a while because there are a lot of parameters...
# but as you can see, these boral models tend to have a ton of parameters

# I've made some custom functions that I think work better for visualling all these parameters
source("code/boral_withOTUs_auxilaryfunctions.R")

plotTrace.lvcoefs(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls) # each line is a parameter associated with an OTU
plotTrace.lvs(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls) # each line is a parameter associated with an OTU
plotTrace.rowcoefs(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls) # each line is a row effect parameter, e.g. there are 3 sites so there are 3 rowID2 lines
plotTrace.Xcoefs(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls) # each line is a parameter associated with an OTU, Xcoef 1 = waterperc, 2 = density, 3 = sizesmall, correspond to column positions in the covariate matrix

# interpreting the traces...
# remember that there's only 1 MCMC chain for all boral models (check out boral help for more details on why)
# your parameter estimate (y) has converged if the estimated value remains stable over the iterations (x) (e.g. chain's slope = 0)

# Calculate convergence statistics -- Geweke convergence diagnostic
# this diagnostic is recommended because of the 1 mcmc chain thing; its discussed on the boral() help page
# Geweke (1992) proposed a convergence diagnostic for Markov chains based on a test for equality of the means of the first and last part of a Markov chain (by default the first 10% and the last 50%). 
# If the samples are drawn from the stationary distribution of the chain, the two means are equal and Geweke's statistic has an asymptotically standard normal distribution
# So, if the p-value is less than 0.05, there is evidence the MCMC chain has not converged

# calculate z
geweke <- geweke.diag(mcmc.obj, frac1=0.1, frac2=0.5)
geweke.z <- geweke[[1]]$z
geweke.df <- data.frame(param = names(geweke.z), z = geweke.z)
# reformat, convert z to pval, do holm pval adjustment for multiple comparisons
geweke.df %>%
  mutate(pval = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
  mutate(adj.pval = p.adjust(pval, method = "holm")) -> geweke.df
# annotate params
geweke.df %>%
  separate(param, into = c("paramType","string"), sep = "\\[", fill = "right") %>%
  separate(string, into = c("OTUnum","string2"), sep = "\\,", fill = "right") %>%
  separate(string2, into = c("paramNum", "drop"), sep = "\\]", fill = "right") %>%
  select(paramType, paramNum, OTUnum, z, pval, adj.pval) -> geweke.df.ann

head(geweke.df.ann) # pval that indicates whether each parameter has converged

geweke.df.ann %>%
  group_by(paramType, paramNum) %>%
  summarize(total = length(adj.pval),
            numConverged = sum(adj.pval > 0.05),
            propConverged = numConverged/total) -> summ
summ # proportion of parameters that have converged by parameter type



# -------------------------------------------------------------------#

# Check residuals

plot(fit.lv)
plot(fit.envlv)


# -------------------------------------------------------------------#

# Visualize as ordinations

par(mfrow=c(2,1))
# unconstrained model
lvsplot(fit.lv, biplot=F)
# constrained model
lvsplot(fit.envlv, biplot=F) # I like to turn off the biplot because the OTU names as vectors get to be overwhelming

# if you'd like to make prettier plots, you can access the scaled lv scores...
lvsplot.vals <- lvsplot(fit.envlv, biplot = TRUE, est = "median", return.vals = TRUE)
lvsplot.vals



# -------------------------------------------------------------------#

# Check out the X coef estimates

#make dfs
OTUId <- row.names(fit.envlv$X.coefs.median)
X.df <- data.frame(OTUId = OTUId, fit.envlv$X.coefs.median)
upper.df <- data.frame(OTUId = OTUId, fit.envlv$hpdintervals$X.coefs[,,"upper"])
lower.df <- data.frame(OTUId = OTUId, fit.envlv$hpdintervals$X.coefs[,,"lower"])
#make long
X.df %>%
  gather(key = "term", value = "X", -OTUId) -> X.df.l
upper.df %>%
  gather(key = "term", value = "X.upper", -OTUId) -> upper.df.l
lower.df %>%
  gather(key = "term", value = "X.lower", -OTUId) -> lower.df.l
#join dfs
X.df.l %>%
  left_join(upper.df.l) %>%
  left_join(lower.df.l) -> Xcoefs.df

head(Xcoefs.df)
# lower and upper refer to the bounds of the 95% highest posterior density (HPD) credible interval


# -------------------------------------------------------------------#

# Check out correlation between OTUs due to...
# (i) similarities in the response to explainatory variables (ie shared enviromental response); see get.enviro.cor() help for more info
# (ii) residual variation accounted for by latent variables; see get.residual.cor() help for more info

enviro.cor <- get.enviro.cor(fit.envlv, prob = 0.95)
residual.cor <- get.residual.cor(fit.envlv, prob = 0.95)

# these objects contain a lot of info, so I wrote a function to extract just the info I want...
source("code/boral_withOTUs_auxilaryfunctions.R")

enviro.cor.df <- extract_cors(corobj = enviro.cor, corType = "enviro")
residual.cor.df <- extract_cors(corobj = residual.cor, corType = "residual")
enviro.cor.df %>%
  left_join(residual.cor.df) -> cor.df

head(cor.df)
dim(cor.df)
# this is a dataframe that holds infor each pairwise OTU combination (rows), i.e. there are 1485 unique pairs of 55 OTUs
# the shared environment and residual correlation values are there for each OTU pair,
# along with an categorical factor that indicates whether the correlation is "significant", ie 95% highest posterior interval does not contain zero

# plot the distribution of correlations in your study

#signif factors
signif.order <- c("no","positive","negative")
cor.df$enviro_signif <- factor(cor.df$enviro_signif, levels = signif.order)
cor.df$residual_signif <- factor(cor.df$residual_signif, levels = signif.order)

#enviro
cor.df %>%
  group_by(enviro_signif) %>%
  summarize(num = length(enviro_signif)) -> tmp
tmp
total <- sum(tmp[,2])
percPos <- round((tmp[tmp$enviro_signif=="positive", "num"]/total)*100, digits = 0)
percNeg <- round((tmp[tmp$enviro_signif=="negative", "num"]/total)*100, digits = 0)
p.env<-ggplot(cor.df, aes(x = enviro_cor, fill = enviro_signif)) + 
  geom_histogram() +
  geom_vline(xintercept = 0, linetype=2) +
  xlab("Shared environment correlation") + ylab("Frequency") +
  annotate("text", y=150, x=-0.5, color = 4, label=paste(percNeg, "%", sep="")) +
  annotate("text", y=150, x=0.5, color = 3, label=paste(percPos ,"%", sep=""))

# Shared enviromental correlation between two OTUs suggests that the OTU abundances either respond in the same way (positively covary) or in opposite ways (negatively covary) to the measured X vars included in the model
# This output doesn't mean too much since the model did not converge, but if it did you could say that... 
# Many OTU pairs either occupy either very similar (positive correlations, 15%) or very different (negative correlations, 14%) wood habitat as characterized by wood size, waterperc, and density
p.env

#residual
cor.df %>%
  group_by(residual_signif) %>%
  summarize(num = length(residual_signif)) -> tmp
total <- sum(tmp[,2])
percPos <- round((tmp[tmp$residual_signif=="positive", "num"]/total)*100, digits = 0)
percNeg <- round((tmp[tmp$residual_signif=="negative", "num"]/total)*100, digits = 0)
p.res<-ggplot(cor.df, aes(x = residual_cor, fill = residual_signif)) + 
  geom_histogram() +
  geom_vline(xintercept = 0, linetype=2) +
  xlab("Residual correlation") + ylab("Frequency") +
  annotate("text", y=150, x=-0.0005, color=4, label=paste(percNeg,"%", sep="")) +
  annotate("text", y=150, x=0.0005, color=3, label=paste(percPos,"%", sep=""))

# Residual correlation between two OTUs suggests that there was some variable(s) not accounted for in the model that influenced the OTU abundances to positively/negatively covary
# This output doesn't mean too much since the model did not converge, but if it did you could say that... 
# Most OTU pairs are not residually correlated after accounting for the measured wood covariates included in the model (100 - (14 + 21) %)
# More OTU pairs have positive (21%) than negative (14%) residual correlations, although the magnitude of those "significant" residual correlations are often weak
p.res


# plot a chord diagram to dig deeper OTU correlations
# here's an example using the shared environment correlations
# you'll probably be interested in annotating with more informative OTU names, but this will do for now...

# filter out the OTU pairs that are not significantly correlated and color based on whether the correlation is positive or negative
cor.df %>%
  filter(enviro_signif != "no") %>%
  select(otu1, otu2, enviro_cor) %>%
  mutate(col = ifelse(enviro_cor > 0, 
                      "black","red")) %>%
  arrange(desc(enviro_cor)) -> env.frame

# start plot
par(mfrow=c(1,1))
circos.par(track.height=.5)
chordDiagram(env.frame[,1:3], 
             col = env.frame$col, 
             annotationTrack = c("grid"), 
             annotationTrackHeight = c(0.01,0.01),
             preAllocateTracks = 1,
             grid.col= 1, reduce = 0)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .001, sector.name, 
              facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5), cex=.7)
}, bg.border = NA)
  
  
  
  
  
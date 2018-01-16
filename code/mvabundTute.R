
### MVABUND package, 2014 updates ###
# David Warton
# School of Mathematics and Statistics, UNSW Australia
# prepared for Eco-Stats lab, September 2014, at the last minute

### What does mvabund do? ###
# Analyses multivariate data (especially abundance of presence-absence data) using simultaneous
# univariate models and design-based inference.
# The main functions are manyglm, which fits a GLM to each response variable, and anova/summary,
# which use row-resampling for valid multivariate inference (i.e. taking into account
# correlation between variables)
# Designed specially for multivariate abundance data in ecology, species-by-site stuff, which
# has two key properties that need to be dealt with:
# (1) strong mean-variance relationship.
# (2) correlation between response variables (e.g. due to species interaction)

### Why is mvabund better than using PRIMER, PC-ORD, etc? ###
# A few reasons:
# - It deals with property (1) as well as (2). Previous methods don't take (1) seriously enough
#   and (at best) try to transform the problem away. Doesn't work when you have lots of zeros.
# - Some important ecological questions imply a model, so you need a model to answer them
#   properly. e.g. the study of composition, interactions, random effects...
# - This all means better statistical properties (good power and valid Type I for more problems)
# - Interpretability... when using models your assumptions are explicitly stated, parameters in
#   your model have explicit meaning. Multivariate stuff is always hard to interpret, we need
#   methods that are more interpretable to make some sense of what is going on.
# - Flexibility... can be generalised to handle all sorts of different problems (different data
#   types, different questions). eg use an offset to account for changes in sampling intensity
# - Diagnostics... really you need to always look at your data and use it to check if the
#   analysis you did was right for your data. This is easier with models like GLMs where your
#   assumptions are clear and many assumption-checking tools are established.
# - and more...

### Why is mavabund worse? ###
# - It is usually slower, although for big datasets (with many observations) it can be faster
#   (it scales linearly with observations not quadratically).
# - In some respects there is less functionality (e.g. no ordination function at the moment).
#   The best course of action here is to harass David and be patient, or help us code!
# - GLMs aren't perfect, there are some traps for young players. Like using summary on pres/abs
#   data, can give weird answers (Wald statistics can be weird). When in doubt, anova is safer.

## Example - Tikus Islands data

library(mvabund)

# subset to the first two years
data(tikus)
is8183 = tikus$x$time=="81"|tikus$x$time=="83" #this tells me which rows are 1981 and 1983
dat = tikus$abund[is8183,]
time = factor(tikus$x$time[is8183]) # a factor telling us the year of sampling
id = factor(tikus$x$rep[is8183]) # a factor telling us the transect that was sampled

# get rid of singletons/doubletons/tripletons/etc (almost never informative!)
nPres = apply(dat>0,2,sum)
datSmall = dat[,nPres>3] # gets rid of anything with less than 3 presences
coverDat = mvabund(datSmall) # so R knows to treat coverDat as multivariate abundances

# plot the data
plot(coverDat~time)
# this function  doesn't work if the window is too small (so resize it if error).
# By default, only the 12 most abundant species are plotted, can be changed (see ?plot.mvabund).
# Note the striking pattern - by 1983, nearly all the coral died!

# Now look at the mean-variance plot
meanvar.plot(coverDat~time,legend=T)
# As the mean increases, so does the variance. Variance ranges over a factor of >1,000
# (most datasets I see have more variation than this). Also points not shown: mean=0, var=0.

# Does transformation fix the problem?
meanvar.plot(log(coverDat+1)~time,legend=T)
# No. You always get the trail off in the bottom left (and always have points at mn=0, vr=0).
# This happens because of all the zeros in the data, and is typical of this type of data.

### "no frills" manyglm analysis ###

# Problem though - manyglm only handles counts (neg binomial, poisson, binomial) not cover.

# Let's start with presence-absence data, which manyglm can handle
presDat = coverDat
presDat[presDat>1] = 1 #turning into pres/abs

# Now look at the mean-variance plot of pres-absence
meanvar.plot(presDat~time,legend=T,log="")
# Note the quadratic up-down thing - pres/absence data ALWAYS looks like this -> it is binomial.

# fit a model assuming main effects of time and transect ID
ft=manyglm(presDat~id+time,family="binomial")
plot(ft)
# no pattern on plot suggests no evidence against model. Try fitting again:
plot(ft)
# notice things changed a bit, because residuals have some randomness in them.

# do an anova, we care mostly about time effect:
an = anova(ft, show.time="all", p.uni="adjusted")
print(an$table)
# Significant effect of time (no kidding)
# If you have a bigger dataset, it takes longer. Can beat this by parallelising and combining
# results, by doing initial analysis with less resamples, eg nBoot=25, or by going to lunch.

# which are the "top 10" species in terms of change with time:
s = sort(an$uni.test[3,],decreasing=T,index.return=T)[1:10]
# an$uni.test stores univariate test stats, 3rd row has change in deviance due to time
s$ix[1:10] # the column numbers of the "top 10" most impacted spp
dimnames(presDat)[[2]][s$ix[1:10]] #the names of the top 10 spp

plot(coverDat~time,var.subset=s$ix[1:10])
# a plot of the top 10. Notice all the circles on the right and trangles on the left - in each
# case, a spp present at most transects has disappeared from all in 1983.

# percentage of change in deviance due to this "top ten":
sum(an$uni.test[3,s$ix[1:10]]) / sum(an$uni.test[3,])
sum(an$uni.test["time",s$ix[1:10]]) / sum(an$uni.test[3,]) #same thing
s$ix[1:10] #column numbers of the top 10 most impacted species
# only about 50% of change in deviance is due to these ten, so we are missing part of the story
# if we focus just on these.
# (The part we are missing is that the other 19 spp mostly disappeared too!)

# You can look at coefficients and their se's too, but in this case, that is pretty meaningless
# as most spp just went to zero.


### the manyany function ###

# manyglm only fits a few common glm's, you can get the rest using the manyany function.
# It can (in principle) take pretty much any univariate model (with predict argument)
# The downside is that it is MUCH slower. Take lunch and then a siesta.

?manyany
# Did this line work? If you get an error, maybe you have an old version of mvabund?
citation("mvabund")
# Does this say you have version 3.9 or more? If not, you'll need to download it!

# in manyany you have to first specify your fitting function ("glm", "gam", "glmer", ...)
# next you specify a matrix or mvabund object with the response in it.
# Then everything else.
# But you need to have a data argument with all your predictors in it.

# e.g. to fit the same model as above but using manyany:
X = data.frame(id=id,time=time)
ftAny = manyany("glm",presDat,presDat~id+time,family="binomial",data=X)
# This gives pretty much the same results as ftAny.
# e.g. anova looks the same, but requires a null object before it will run.
ftNullAny = manyany("glm",presDat,presDat~id,family="binomial",data=X)
anova(ftNullAny,ftAny)
# Note that this produces pretty much the same results as before, but takes heaps longer!
# Warning: default nBoot is now 99 not 999, because things are so much slower.
# for publication you should ramp this up to 999, then go to lunch.


### cool stuff you can do with the manyany function ###

# COMPOSITONAL models: you can test for changes in relative not absolute aboundance.
# Other packages (PRIMER etc) are not good at doing this.

ftAny = manyany("glm", presDat, presDat~id+time, family="binomial", data=X, composition=TRUE)
ftNullAny = manyany("glm", presDat, presDat~id, family="binomial", data=X, composition=TRUE)
anova(ftNullAny,ftAny,nBoot=9)
# you need a lot more than 9 resamples for a proper test, but you'd need to go for coffee!
# I have thoughts on how to make this faster but haven't tried them out yet.

# CHANGE LINK FUNCTION: you can use anything from glm not yet in manyglm.
ftProbit = manyany("glm", presDat, presDat~id+time, family=binomial("probit"),
                   data=X,composition=TRUE)

# CHANGE FAMILY: for cover or biomass data, you could try using the Tweedie distribution:
library(statmod)
library(tweedie)
ftTw = manyany("glm", coverDat, y~id+time, family=tweedie(var.power=1.5, link.power=0),
               var.power=1.5, data=X, composition=TRUE)
# Family functions need to be added to code by David on a case-by-case basis but easy to do

# FIT A GAM: for example, using the spider data. But note that this is "data-hungry" so not
#  worth much if you have a small number of sites or for spp with heaps of zerps

# load a different dataset with continuous predictors:
data(spider)
presSpider <- spider$abund
presSpider[presSpider>1] = 1
X <- data.frame(spider$x)
require(mgcv)
require(MASS)
ftGAM = manyany("gam", presSpider, y~s(soil.dry), family=binomial(link="cloglog"), data=X)
ftGLM = manyany("glm", presSpider, y~poly(soil.dry,2), family=binomial(link="cloglog"), data=X)
print( sum(BIC(ftGAM)) )
print( sum(BIC(ftGLM)) )
# The GAM has a heaps smaller BIC than the quadratic model so is a better fit to the data

# MIXED EFFECTS MODELS: can be fitted via lme4, but no testing via anova function as yet.

# BLOCK RESAMPLING: can use the argument block=... to specify a block factor used to resample
#  blocks of rows instead of resampling individual rows. This can be used for repeated measures,
#  or to fit trait models (but requires a bit of data manipulatation first, another time?) 


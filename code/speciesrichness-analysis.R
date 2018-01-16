library(boral) ## install.packages("boral") if not done previously
source("speciesrichness-auxilaryfunctions.R") ## R script containing functions for calculating species richness

## The tree dataset is not available publicly. In the below code however, it is assumed that the data has been loaded and consists of four objects: 1) a presence-absence response matrix called Y for training data, 2) a matrix of two covariates called X with column names DDEG (growing degree days) and PRSU (summer preciptation) for training data, 3) a presence-absence response matrix called Ytest for test data, 2) a matrix of two covariates called Xtest for test data. 
## To affirm, the full dataset has n=6118 sites and p=51 tree species

## Create quadratic terms for each covariate for training data
X2 <- cbind(poly(X$DDEG,2), poly(X$PRSU,2)) 

## Create quadratic terms for each covariate for the test data
Xtest2 <- cbind(predict.poly(poly(X$DDEG,2), Xtest$DDEG), predict.poly(poly(X$PRSU,2), Xtest$PRSU))
Xtest2 <- Xtest2


##########
## Fit LVM with 2 latent variables to training data, and predict spp richness on to the test sites
## This arguments below can be altered any number of LVs by changing the num.lv argument, including num.lv = 0 for independence model
##########
fit.lv2 <- boral(y = Y, X = X2, family = "binomial", trial.size = 1, num.lv = 2, save.model = TRUE, n.iteration = 60000, n.thin = 50, n.burnin = 10000, calc.ics = F, hypparams = c(20,20,20,20))

fit.mcmcBase <- 
all.mcmc <- mcmc(fit.lv2$jags.model$BUGSoutput$sims.matrix, start = 1, thin = fit.lv2$jags.model$BUGSoutput$n.thin) 
all.lv2.rich <- matrix(NA,nrow(all.mcmc),nrow(Xtest2))
for(t in 1:nrow(all.mcmc)) {
	cw.fit <- make.fit(cw.mcmc.samples = all.mcmc[t,], p = ncol(Ytest), n = nrow(Y))
	all.lv2.rich[t,] <- calc.indices(family = "binomial", p = ncol(Ytest), num.lv = fit.lv2$num.lv, cw.fit = cw.fit, test.y = Ytest, newX = Xtest2, B = 1000)[2,]
	}

	
## Calculate 95% credible intervals for species richness
ci.lv2.rich <- apply(all.lv2.rich,2,quantile,prob = c(0.025,0.975))	


## Calculate true species richess for test sites
true.rich <- apply(Ytest>0,1,sum)


## Determine how many of the credible intervals contain the true species richness
inout.fun <- function(t, ci, true) { out <- ifelse(ci[1,t] < true[t] & ci[2,t] > true[t], 1, 0) }
table(sapply(1:ncol(ci.nolv.rich), inout.fun, ci = ci.lv2.rich, true = true.rich))


## Summary of the interval widths
lv2.ciwidth <- summary(ci.lv2.rich[2,]-ci.lv2.rich[1,])


##########
## Now take the LVM fitted above with 2 latent variables, and calculate conditional species richness onto the training sites
##########
all.mcmc <- mcmc(fit.lv2$jags.model$BUGSoutput$sims.matrix, start = 1, thin = fit.lv2$jags.model$BUGSoutput$n.thin) 
all.lv2.rich <- array(NA,dim=c(nrow(all.mcmc),3,nrow(X2)))
for(t in 1:nrow(all.mcmc)) {
	cw.fit <- make.fit(cw.mcmc.samples = all.mcmc[t,], p = ncol(Y), n = nrow(Y))
	all.lv2.rich[t,,] <- calc.cond.indices(family = "binomial", p = ncol(Ytest), num.lv = 2, cw.fit = cw.fit, test.y = Y, test.X = X2, B = 1000)
	}

fit.lv2$jags.model <- NULL	


## Calculate true species richess 
true.rich <- apply(Ytest>0,1,sum)


## Average over the MCMC samples. Order: True, conditional, marginal
lv2.rich.mean <- rbind(apply(all.lv2.rich[,1,],2,mean), apply(all.lv2.rich[,2,],2,mean),apply(all.lv2.rich[,3,],2,mean))


## Mean squared error for species richness (should be smaller for conditional values)
## Scaled by multiplying by ten for easier reading
as.matrix(dist(lv2.rich.mean))[,1]/length(true.rich)*10 


## Calculate 95% credible intervals and their summary for conditional species richness
ci.lv2.cond.rich <- apply(all.lv2.rich[,2,],2,quantile,prob = c(0.025,0.975))
lv2.ciwidth <- summary(ci.lv2.cond.rich[2,]-ci.lv2.cond.rich[1,])


## Calculate 95% credible intervals and their summary for marginal species richness (should be wider than conditional values)
ci.lv2.marg.rich <- apply(all.lv2.rich[,3,],2,quantile,prob = c(0.025,0.975)) 
lv2.ciwidth <- summary(ci.lv2.marg.rich[2,]-ci.lv2.marg.rich[1,])



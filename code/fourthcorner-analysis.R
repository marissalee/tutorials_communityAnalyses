rm(list = ls())
library(ade4)
source("./scripts/fourthcorner-auxilaryfunctions.R") ## R script containing functions for fitting models for fourth-corner analysis

data(aravo)
str(aravo)
#addresses the potential impact of snowmelt dates on plant communities
#alpine plants at 75 sites in Aravo, SE France
#use presence/absence to predict the probability of an observed presence of species j in site i as a quadratic function of mean snowmelt date, SLA, and 2 latent vars
# prob(mij) ~ ai + b0j + snowi(b1) + snow^2i(b2) + snowi(SLAj)(b3) + snowi(bj) + zi(lambdaj)
# ai and bj are random effects

## Convert to presence-absence response
aravo.pa <- (aravo$spe>0) 
aravo.pa[aravo.pa == FALSE] = 0; 
aravo.pa[aravo.pa == TRUE] = 1

## Consider only species with more than four presences
sel.spp <- which(colSums(aravo.pa) > 4) 
aravo.pa <- aravo.pa[,sel.spp]


##########
## Setup some stuff that will be needed for all the models to be fitted in JAGS
##########
gen.inits <- function() {  
	Tau <- rWishart(1,ncol(aravo.pa)+1,diag(ncol(aravo.pa)))[,,1]
	Sigma <- solve(Tau)
	Z <- abs(t(rmvnorm(nrow(aravo.pa),rep(0,ncol(aravo.pa)),Sigma)))
	Z <- ifelse(as.matrix(aravo.pa), Z, -1 * Z)
	list("Z" = Z) 
	}
	

##########
## Fit independence model with no latent variables but includes interaction term
##########	
basic.nolvmodfile <- make.boralmodel.fourthcorner(family = "binomial", 
                                                  X.eff = "ran", site.eff = "ran", 
                                                  num.lv = 0, fc.eff = T, 
                                                  filename = "nolvint.txt")
jags.data <- list(y = aravo.pa, n = nrow(aravo.pa), 
                  p = ncol(aravo.pa), traits = as.vector(scale(aravo$traits$SLA[sel.spp])), 
                  enviro = as.vector(scale(aravo$env$Snow)))
fit.nolvint <- jags(data = jags.data, 
                    inits = gen.inits, 
                    parameters.to.save = c("mean.X.eff", "X.params", 
                                           "quad.X.params", "sigma2.X.eff", 
                                           "sigma2.site", "all.params",
                                           "inter.params","site.params"), 
                    model.file = "nolvint.txt", n.chains = 1, n.iter = 55000, n.burnin = 5000, n.thin = 50)

## Look at 95% credible interval for the interaction term
intervals.nolv <- HPDinterval(as.mcmc(fit.nolvint))[[1]]
round(intervals.nolv["inter.params",],3)
#yes, the interaction term differs from 0 -- there is an interaction between environment and traits that drive abundance observations

##########
## Fit LVM with two latent variables and includes interaction term
##########
basic.modfile <- make.boralmodel.fourthcorner(family = "binomial", 
                                              X.eff = "ran", 
                                              site.eff = "ran", 
                                              num.lv = 2, 
                                              fc.eff = T, filename = "lvint.txt")
jags.data <- list(y = aravo.pa, 
                  n = nrow(aravo.pa), 
                  p = ncol(aravo.pa), 
                  traits = as.vector(scale(aravo$traits$SLA[sel.spp])), 
                  enviro = as.vector(scale(aravo$env$Snow)), num.lv = 2)
fit.lvint <- jags(data = jags.data, 
                  inits = gen.inits, 
                  parameters.to.save = c("mean.X.eff","quad.X.params",
                                         "X.params", "sigma2.X.eff", 
                                         "sigma2.site","all.params",
                                         "inter.params","lvs"), 
                  model.file = "lvint.txt", n.chains = 1, 
                  n.iter = 55000, n.burnin = 5000, n.thin = 50)


## Look at 95% credible interval for the interaction term -- this should be wider than the interval obtained from the independence model above
intervals.lv <- HPDinterval(as.mcmc(fit.lvint))[[1]]
round(intervals.lv["inter.params",],3)


## Extract posterior means of latent variables and their coefficients from MCMC output
all.mcmc <- mcmc(fit.lvint$BUGSoutput$sims.matrix, start = 1, thin = fit.lvint$BUGSoutput$n.thin) 
params.mean <- colMeans(all.mcmc)
fit.lvint$lv.mean <- matrix(params.mean[grep("lvs",colnames(all.mcmc))],nrow(aravo.pa))
fit.lvint$lv.coefs.mean <- matrix(params.mean[grep("all.params",colnames(all.mcmc))],ncol(aravo.pa))
rm(all.mcmc)


##########
## Fit a pure LVM including only two latent variables 
##########
basic.modfile <- make.boralmodel.fourthcorner(family = "binomial", 
                                              X.eff = "none", 
                                              site.eff = "ran", 
                                              num.lv = 2, 
                                              fc.eff = F, 
                                              filename = "purelv.txt")
jags.data <- list(y = aravo.pa, n = nrow(aravo.pa), 
                  p = ncol(aravo.pa), num.lv = 2)
fit.purelv2 <- jags(data = jags.data, 
                    inits = gen.inits, 
                    parameters.to.save = c("site.params","sigma2.site", 
                                           "all.params","lvs"), 
                    model.file = "purelv.txt", n.chains = 1, 
                    n.iter = 55000, n.burnin = 5000, n.thin = 50)

## Extract posterior means of latent variables and their coefficients from MCMC output
all.mcmc <- mcmc(fit.purelv2$BUGSoutput$sims.matrix, start = 1, thin = fit.purelv2$BUGSoutput$n.thin) 
params.mean <- colMeans(all.mcmc)
fit.purelv2$lv.mean <- matrix(params.mean[grep("lvs",colnames(all.mcmc))],nrow(aravo.pa))
fit.purelv2$lv.coefs.mean <- matrix(params.mean[grep("all.params",colnames(all.mcmc))],ncol(aravo.pa))
rm(all.mcmc)


## Compare the trace of Lambda%*%t(Lambda) for the pure LVM vs the fourth corner LVM -- it should drop dramatically reflecting the amount of covariation explained by the fourth corner model
purelv2.rescov <- (fit.purelv2$lv.coefs.mean[,2:3]%*%t(fit.purelv2$lv.coefs.mean[,2:3])) #pure LVM
lv2.rescov <- (1*fit.lvint$lv.coefs.mean[,2:3]%*%t(1*fit.lvint$lv.coefs.mean[,2:3])) #4thcorner LVM

sum(diag(purelv2.rescov))
sum(diag(lv2.rescov))

explain.rat <- (sum(diag(purelv2.rescov))-sum(diag(lv2.rescov)))/sum(diag(purelv2.rescov))
explain.rat #18% drop in the trace ... how to interpret this?

##########
## For the model with two latent variables and interaction term (fit.lvint), construct a residual correlation plot
########## 

## Function to calculate the residual correlation matrix from a LVM
get.residual.corv2 <- function(fit.boral, adjust = 1) {
	fit.mcmcBase <- fit.boral$BUGSoutput
	mcmc.runs <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = fit.mcmcBase$n.thin) 
	rm(fit.mcmcBase)

	y <- fit.boral$y
	n.species <- ncol(y)
	n.sites <- nrow(y)
	
	num.lv <- fit.boral$num.lv
	Tau.arr <- matrix(NA,nrow(mcmc.runs),n.species^2)
	Tau.cor.arr <- matrix(NA,nrow(mcmc.runs),n.species^2)
	all.sigmas2 <- matrix(NA,nrow(mcmc.runs),n.species)
	for(t in 1:nrow(mcmc.runs)) { 
		lvs <- matrix(mcmc.runs[t,grep("lvs",colnames(mcmc.runs))],n.sites,byrow=F)
		lv.coefs <- matrix(mcmc.runs[t,grep("all.params",colnames(mcmc.runs))],n.species,byrow=F)
		lv.coefs2 <- lv.coefs[,2:(num.lv+1)]
		Tau.mat <- lv.coefs2%*%t(lv.coefs2) + diag(n.species)
		Tau.arr[t,] <- as.vector(Tau.mat) 
		
		Tau.cor.mat <- cov2cor(Tau.mat)
		Tau.cor.arr[t,] <- as.vector(Tau.cor.mat) 
		}
		
	## Average/Median over the MCMC samples
	Tau.mat.mean <- sig.Tau.mat.mean <- matrix(apply(Tau.arr,2,mean),n.species,byrow=F)
	Tau.mat.median <- sig.Tau.mat.median <- matrix(apply(Tau.arr,2,median),n.species,byrow=F)
	Tau.cor.mean <- sig.Tau.cor.mean <- matrix(apply(Tau.cor.arr,2,mean),n.species,byrow=F)
	Tau.cor.median <- sig.Tau.cor.median <- matrix(apply(Tau.cor.arr,2,median),n.species,byrow=F)
		
	get.cor.intervals <- HPDinterval(as.mcmc(Tau.arr), prob = 0.95)	
	id.sign.cors <- which(get.cor.intervals[,1] > 0 | get.cor.intervals[,2] < 0)	
	sig.Tau.mat.mean[-id.sign.cors] <- 0
	sig.Tau.mat.median[-id.sign.cors] <- 0
	sig.Tau.cor.median[-id.sign.cors] <- 0	
	sig.Tau.cor.mean[-id.sign.cors] <- 0	
		
	return(list(cov.mean = Tau.mat.mean, cov.median = Tau.mat.median, sig.cor.mean = sig.Tau.cor.mean, sig.cor.median = sig.Tau.cor.median, cor.mean = Tau.cor.mean, cor.median = Tau.cor.median))
	}	


fit.lvint$y <- aravo.pa; 
fit.lvint$num.lv <- 2

lv2.cor <- get.residual.corv2(fit.lvint)


## Plot of residual correlation matrix (posterior mean estimator, and reordered by first principal component)
library(corrplot)
lv2.cor$reorder.cor.mean <- corrMatOrder(lv2.cor$cor.mean, order = "FPC", hclust.method = "average")
rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- sapply(strsplit(aravo$spe.names[sel.spp]," "),function(x) paste(x[1:2],collapse=" "))

par(cex = 1, cex.main = 1.5)
corrplot(lv2.cor$cor.mean[lv2.cor$reorder.cor.mean,lv2.cor$reorder.cor.mean], diag = F, type = "lower", title = "Residual Correlation Matrix from LVM", mar = c(1,1,3,1), method = "color", tl.srt = 45, tl.cex = 0.5)



##########
## Unconstrained ordination biplot -- results can vary depending on hypparameters
##########
## Change a couple of names to make it consistent with correlation plot above.
colnames(aravo.pa)[which(colnames(aravo.pa) == "Card.alpi")] <- "Card.bell"
colnames(aravo.pa)[which(colnames(aravo.pa) == "Cera.stri")] <- "Card.arve"


## Extract posterior means of latent variables and their coefficients from MCMC output
all.mcmc <- mcmc(fit.purelv2$BUGSoutput$sims.matrix, start = 1, thin = fit.purelv2$BUGSoutput$n.thin) 
params.mean <- colMeans(all.mcmc)
fit.purelv2$lv.mean <- matrix(params.mean[grep("lvs",colnames(all.mcmc))],nrow(aravo.pa))
fit.purelv2$lv.coefs.mean <- matrix(params.mean[grep("all.params",colnames(all.mcmc))],ncol(aravo.pa))
rm(all.mcmc)
rownames(fit.purelv2$lv.coefs.mean) <- colnames(aravo.pa)


## Adjust scaling so that both the site and species ordinations are on the same scale. 
fit.purelv2$lv.mean2 <- fit.purelv2$lv.mean*matrix(sqrt(colSums(fit.purelv2$lv.coefs.mean[,2:3]^2))/sqrt(colSums(fit.purelv2$lv.mean^2)),nrow(aravo.pa),2,byrow=T) 
fit.purelv2$lv.coefs.mean2 <- fit.purelv2$lv.coefs.mean[,2:3]*matrix(sqrt(colSums(fit.purelv2$lv.mean2^2))/sqrt(colSums(fit.purelv2$lv.coefs.mean[,2:3]^2)),ncol(aravo.pa),2,byrow=T) 


library(RColorBrewer)
snow.cols <- brewer.pal(length(unique(aravo$env$Snow)),"PuBuGn")

plot(fit.purelv2$lv.mean2, type = "n", xlab = "LV1", ylab = "LV2", main = "(a) Unconstrained ordination", xlim = c(-3,3), ylim = c(-3,3)) ## xlim and ylim may need to be adjusted depending on values obtained from LVM
text(fit.purelv2$lv.mean2, label = aravo$env$Snow/10, col = snow.cols[as.numeric(factor(aravo$env$Snow))])
## Select the 20 species with the largest L2 norms 
largest.lnorm <- order(apply((fit.purelv2$lv.coefs.mean2)^2,1,sum),decreasing = T)[1:20] 
points(fit.purelv2$lv.coefs.mean2[largest.lnorm,], pch=15, cex = 0.75)
text(fit.purelv2$lv.coefs.mean2[largest.lnorm,1]+ 0.1*sign(fit.purelv2$lv.coefs.mean2)[largest.lnorm,1], fit.purelv2$lv.coefs.mean2[largest.lnorm,2] + 0.1*sign(fit.purelv2$lv.coefs.mean2)[largest.lnorm,2], label = colnames(aravo.pa)[largest.lnorm], cex = 0.75) ## Shift spp names slightly to avoid clash with points
legend.col(col = snow.cols, lev = aravo$env$Snow)


##########
## Residual ordination biplot -- results can vary depending on hypparameters
## A procrustes rotation could be applied so you get the same orientation as the Unconstrained biplot
##########
## Extract posterior means of latent variables and their coefficients from MCMC output
all.mcmc <- mcmc(fit.lvint$BUGSoutput$sims.matrix, start = 1, thin = fit.lvint$BUGSoutput$n.thin) 
params.mean <- colMeans(all.mcmc)
fit.lvint$lv.mean <- matrix(params.mean[grep("lvs",colnames(all.mcmc))],nrow(aravo.pa))
fit.lvint$lv.coefs.mean <- matrix(params.mean[grep("all.params",colnames(all.mcmc))],ncol(aravo.pa))
rm(all.mcmc)

rownames(fit.purelv2$lv.coefs.mean) <- rownames(fit.lvint$lv.coefs.mean) <- sapply(strsplit(aravo$spe.names[sel.spp]," "),function(x) paste(x[1:2],collapse=" "))


## Adjust scaling so that both the site and species ordinations are on the same scale. 
fit.lvint$lv.mean2 <- fit.lvint$lv.mean*matrix(sqrt(colSums(fit.lvint$lv.coefs.mean[,2:3]^2))/sqrt(colSums(fit.lvint$lv.mean^2)),nrow(aravo.pa),2,byrow=T) 
fit.lvint$lv.coefs.mean2 <- fit.lvint$lv.coefs.mean[,2:3]*matrix(sqrt(colSums(fit.lvint$lv.mean2^2))/sqrt(colSums(fit.lvint$lv.coefs.mean[,2:3]^2)),ncol(aravo.pa),2,byrow=T) 

snow.cols <- brewer.pal(length(unique(aravo$env$Snow)),"PuBuGn")

#fit.lvint$lv.mean2 <- procrustes(fit.purelv2$lv.mean2, fit.lvint$lv.mean2, scale=F)$Yrot
#fit.lvint$lv.coefs.mean2 <- procrustes(fit.purelv2$lv.coefs.mean2, fit.lvint$lv.coefs.mean2, scale=F)$Yrot
plot(fit.lvint$lv.mean2, type = "n", xlab = "LV1", ylab = "LV2", main = "(b) Residual ordination", ylim = c(-3,3), xlim = c(-3,3)) ## xlim and ylim may need to be adjusted depending on values obtained from LVM
text(fit.lvint$lv.mean2, label = aravo$env$Snow/10, col = snow.cols[as.numeric(factor(aravo$env$Snow))])
## Select the 20 species with the largest L2 norms 
largest.lnorm <- order(apply((fit.lvint$lv.coefs.mean2)^2,1,sum),decreasing = T)[1:20]
points(fit.lvint$lv.coefs.mean2[largest.lnorm,],pch=15, cex = 0.75)
text(fit.lvint$lv.coefs.mean2[largest.lnorm,1]+ 0.1*sign(fit.lvint$lv.coefs.mean2)[largest.lnorm,1], fit.lvint$lv.coefs.mean2[largest.lnorm,2] + 0.1*sign(fit.lvint$lv.coefs.mean2)[largest.lnorm,2], label = colnames(aravo.pa)[largest.lnorm], cex = 0.75) ## Shift spp names slightly to avoid clash with points
legend.col(col = snow.cols, lev = aravo$env$Snow)





#############
## UNUSED STUFF
#############
# ## Perform a Likelihood Ratio test using Bayesian MCMC samples (admittedly a bit dodgy, and prone sometimes to sampling error) -- anyway, the test ##########
# statistic for independence models (lrt.nolv) should be much larger compared to the LVMs (lrt.lv)
# fit.nolvint$logL <- calc.logL(jagsfit = fit.nolvint, num.lv = 0, y = aravo.pa, X = X2, inter.array = inter.array)
# fit.nolvnoint$logL <- calc.logL(jagsfit = fit.nolvnoint, num.lv = 0, y = aravo.pa, X = X2, inter.array = inter.array)
# 
# fit.lvint$logL <- calc.logL(jagsfit = fit.lvint, num.lv = 2, y = aravo.pa, X = X2, inter.array = inter.array)
# fit.lvnoint$logL <- calc.logL(jagsfit = fit.lvnoint, num.lv = 2, y = aravo.pa, X = X2, inter.array = inter.array)
# 
# lrt.nolv <- -2*(fit.nolvnoint$logL$max.logL-fit.nolvint$logL$max.logL)
# lrt.lv <- -2*(fit.lvnoint$logL$max.logL-fit.lvint$logL$max.logL)
# 

# ##########
# ## Fit independence model with no latent variables and interaction term
# ##########
# basic.nolvmodfile <- make.jagsboralnullmodel.fourthcorner.noint(n = nrow(aravo.pa), p = ncol(aravo.pa))
# fit.nolvnoint <- jags(data = jags.data, inits = gen.inits, parameters.to.save = c("X.params","all.params"), model.file = "jagsboralnullmodel.txt", n.chains = 1, n.iter = 120000, n.burnin = 20000, n.thin = 100)
# 
# 
# ##########
# ## Fit LVM with two latent variables but excludes interaction term
# ##########
# basic.modfile <- make.jagsboralmodel.fourthcorner.noint(n = nrow(aravo.pa), p = ncol(aravo.pa))
# fit.lvnoint <- jags(data = jags.data, inits = gen.inits, parameters.to.save = c("X.params","all.params","lvs"), model.file = "jagsboralmodel.txt", n.chains = 1, n.iter = 120000, n.burnin = 20000, n.thin = 100)

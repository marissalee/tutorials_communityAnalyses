##############
## Auxilary functions for calculating species richness
##############
library(mvtnorm)

## Predicts onto test sites, using the marginal log-likelihood i.e., averaging of the latent variable distribution. 
calc.indices <- function(family, p, num.lv, cw.fit, test.y, newX = NULL, B = 1000) {
	n <- nrow(newX); 
	if(missing(newX)) stop("newX is required")
	
	true.spp.rich <- apply(test.y>0,1,sum)
	sim.y <- matrix(NA,nrow(test.y),ncol(test.y))
	if(num.lv == 0) {
		cw.eta <- matrix(1,n,1)%*%t(cw.fit$lv.coefs[,1]) + as.matrix(newX)%*%t(cw.fit$X.coefs) ## Calculate predicted mean response 
		for(j in 1:p) { sim.y[,j] <- rbinom(n,1,p=pnorm(cw.eta[,j])) } ## Simulate a predicted response matrix
		est.spp.rich <- apply(sim.y>0,1,sum) ## Predict species richness based on response matrix
		}

	if(num.lv > 0) {
		est.spp.rich.all <- matrix(n,B)
		for(t in 1:B) {
			sub.mc.lvs <- rmvnorm(n,rep(0,num.lv)) ## Randomly sample a set of latent variables from a large number generated from N(0,1) distributions
			cw.eta <- cbind(1,sub.mc.lvs,newX)%*%t(cbind(cw.fit$lv.coefs[,1:(1+num.lv)],cw.fit$X.coefs)) ## Calculate predicted mean response 
			for(j in 1:p) { sim.y[,j] <- rbinom(n,1,p=pnorm(cw.eta[,j])) } ## Simulate a predicted response matrix
			est.spp.rich.all[,t] <- rowSums(sim.y>0) } ## Predict species richness based on response matrix
	
		est.spp.rich <- apply(est.spp.rich.all,1,median) ## Use median instead of mean to avoid non-integer predictions
		}
	
	spp.richness <- rbind(true.spp.rich,est.spp.rich); 
	rownames(spp.richness) <- c("true","estimated")
	
	return(spp.richness)
	}


## Predicts onto training sites used to fit the model i.e., conditional on the latent variables observed
calc.cond.indices <- function(family, p, num.lv, cw.fit, test.y, test.X = NULL, B = 1000) {
	n <- nrow(test.X); 
	if(missing(test.X)) stop("test.X is required")
	
	true.spp.rich <- apply(test.y>0,1,sum)
	if(num.lv == 0) {
		cw.eta <- matrix(1,n,1)%*%t(cw.fit$lv.coefs[,1]) + as.matrix(test.X)%*%t(cw.fit$X.coefs) ## Calculate fitted mean response 
		logL0.fit.nolv <- dbinom(1,1,p=pnorm(cw.eta)) ## Calculate predicted probability of presence for training sites
		est.spp.rich <- rowSums(logL0.fit.nolv) } ## Calculate species richness for training sites

	if(num.lv > 0) {
		cw.eta <- matrix(1,n,1)%*%t(cw.fit$lv.coefs[,1]) + as.matrix(test.X)%*%t(cw.fit$X.coefs) + as.matrix(cw.fit$lvs)%*%t(cw.fit$lv.coefs[,2:(num.lv+1)]) ## Calculate fitted mean response conditional on latent variables
		logL0.fit.nolv <- dbinom(1,1,p=pnorm(cw.eta)) ## Calculate predicted probability of presence for training sites
		est.cond.spp.rich <- rowSums(logL0.fit.nolv) ## Calculate species richness for training sites
			
		est.marg.spp.rich.all <- matrix(n,B)
		for(t in 1:B) {
			sub.mc.lvs <- rmvnorm(n,rep(0,num.lv)) ## Randomly sample a set of latent variables from a large number generated from N(0,1) distributions
			cw.eta <- cbind(1,sub.mc.lvs,newX)%*%t(cbind(cw.fit$lv.coefs[,1:(1+num.lv)],cw.fit$X.coefs)) ## Calculate predicted mean response 
			logL0.fit.nolv <- dbinom(1,1,p=pnorm(cw.eta)) ## Calculate predicted probability of presence for training sites
			est.marg.spp.rich.all[,t] <- rowSums(logL0.fit.nolv) } ## Calculate species richness for training sites
		
		est.marg.spp.rich <- apply(est.spp.rich.all,1,median) ## Use median instead of mean to avoid non-integer predictions
		}
	
	if(num.lv == 0) { spp.richness <- rbind(true.spp.rich,est.spp.rich); rownames(spp.richness) <- c("true","estimated") }
	if(num.lv > 0) { spp.richness <- rbind(true.spp.rich,est.cond.spp.rich,est.marg.spp.rich); rownames(spp.richness) <- c("true","conditional","marginal") } 
	
	return(spp.richness)
	}

	
## Construct list containing key parts of a LVM or independence model	
make.fit <- function(cw.mcmc.samples, p, n = NULL) {
	get.betas <- matrix(cw.mcmc.samples[grep("X.params",names(cw.mcmc.samples))], nrow = p, byrow = F)
	get.lv.coefs <- matrix(cw.mcmc.samples[grep("all.params",names(cw.mcmc.samples))], nrow = p, byrow = F)
	get.lvs <- matrix(cw.mcmc.samples[grep("lvs",names(cw.mcmc.samples))], nrow = n, byrow = F)

	out <- list(lv.coefs = get.lv.coefs, X.coefs = get.betas, lvs = get.lvs)
	return(out)
	}
	

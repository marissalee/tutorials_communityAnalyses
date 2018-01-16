############
## Auxilary functions for fourth corner analysis for the Aravo dataset
##############
library(mvtnorm)
library(R2jags)
library(mvabund)
library(vegan); 


calc.logL <- function(jagsfit, num.lv = 0, y, X = NULL, inter.array, B = 1000) {
	n <- nrow(y)
	p <- ncol(y)
	if(num.lv > 0) mcmc.lvs <- rmvnorm(B,rep(0,num.lv))
	
	fit.mcmcBase <- jagsfit$BUGSoutput
	all.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = jagsfit$BUGSoutput$n.thin) 
	out.logL <- numeric(nrow(all.mcmc))
	
	for(t in 1:nrow(all.mcmc)) {
		all.params <- matrix(all.mcmc[t,grep("all.params",colnames(all.mcmc))],p)
		if(!is.null(X)) { X.params <- all.mcmc[t,grep("X.params",colnames(all.mcmc))] }
		inter.params <- all.mcmc[t,grep("inter.params",colnames(all.mcmc))]
	
		etas <- matrix(all.params[,1],n,p,byrow=T)
		if(!is.null(X)) etas <- etas + matrix(X%*%X.params,n,p,byrow=F) 
		if(length(inter.params) != 0) etas <- etas + inter.array[,,1]*inter.params
		
		if(num.lv == 0) { 
			logL <- dbinom(y, 1, prob = pnorm(etas), log = T) 
			out.logL[t] <- sum(logL[is.finite(logL)]) }
			
		if(num.lv > 0) {
			if(t%%1 == 0) cat("Up to iteration",t,"\n")
			for(i in 1:n) { 
				etas.i <- etas[rep(i,B),] + mcmc.lvs%*%t(all.params[,-1])
				logL <- dbinom(y[rep(i,B),], 1, prob = pnorm(etas.i))
				logL[!is.finite(logL)] <- NA
				out.logL[t] <- out.logL[t] + log(mean(apply(logL,1,prod,na.rm=T))) }
			}
		}	
			
	return(list(logL.all = out.logL, max.logL = max(out.logL))) 
	}
		

## Function to write JAGS code for fitting LVM fourth corner models
## Note the tigther priors on LV coefficients (spp scores), to help facilitate convergence for the rare species 
make.boralmodel.fourthcorner <- function(family = "binomial", num.lv = 2, X.eff = "ran", site.eff = "ran", fc.eff = FALSE, hypparams = c(100,20), filename = "jagsboralmodel.txt") {
	mod.general.lv <- paste("model { \n\t ## Likelihood \n\t for(i in 1:n) { for(j in 1:p) {", sep = "")

	if(site.eff != "none") site.effect.string <- paste(" + site.params[i]", sep = "")
	if(site.eff == "none") site.effect.string <- NULL

	if(fc.eff) fourthcorner.effect <- paste(" + inter.params*enviro[i]*traits[j]", sep = "")
	if(!fc.eff) fourthcorner.effect <- NULL
		
	if(X.eff == "none") X.effect.string <- NULL
	if(X.eff == "fixed") X.effect.string <- paste(" + enviro[i]*X.params[1] + pow(enviro[i],2)*quad.X.params", sep = "")
	if(X.eff == "ran") X.effect.string <- paste(" + enviro[i]*X.params[j] + pow(enviro[i],2)*quad.X.params", sep = "")

	if(num.lv > 0) lv.effect.string <- paste("+ inprod(all.params[j,2:(num.lv+1)],lvs[i,])", sep = "")
	if(num.lv == 0) lv.effect.string <- NULL

	mod.general.lv <- c(mod.general.lv, paste("\t\t eta[i,j] <- all.params[j,1] ", site.effect.string, X.effect.string, fourthcorner.effect, lv.effect.string, sep = ""))
			
	if(family == "binomial") {			
		mod.general.lv <- c(mod.general.lv, paste("\t\t Z[i,j] ~ dnorm(eta[i,j],1)", sep = ""))
		mod.general.lv <- c(mod.general.lv, paste("\t\t y[i,j] ~ dbern(step(Z[i,j])) \n \t\t } } \n", sep = "")) }
	
		
	if(num.lv > 0) 
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } } ## Prior distributions for LVs \n", sep = ""))

	mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:p) { all.params[i,1] ~ dnorm(0,", 1/hypparams[1], ") } ## Species intercept \n", sep = ""))
	if(X.eff != "none") 
		mod.general.lv <- c(mod.general.lv, paste("\t quad.X.params ~ dnorm(0,", 1/hypparams[1], ") \n", sep = ""))
	if(X.eff == "fixed") 
		mod.general.lv <- c(mod.general.lv, paste("\t X.params ~ dnorm(0,", 1/hypparams[1], ") } \n", sep = ""))
	if(X.eff == "ran") {
		mod.general.lv <- c(mod.general.lv, paste("\t sigma2.X.eff ~ dgamma(",1/hypparams[1], ",",1/hypparams[1], ") \n\t mean.X.eff ~ dnorm(0,", 1/hypparams[1], ")", sep = "")) 
		mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j] ~ dnorm(mean.X.eff,1/sigma2.X.eff) } \n", sep = "")) }
   
	if(num.lv == 2) {
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { all.params[i,j] <- 0 } } ## Constraints to 0 on upper diagonal", sep = ""))
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,", hypparams[2], ") } ## Sign constraints on diagonal elements", sep = ""))
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 2:num.lv) { for(j in 2:i) { all.params[i,j] ~ dnorm(0,", 1/hypparams[2], ") } } ## Free lower diagonals", sep = ""))
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in (num.lv+1):p) { for(j in 2:(num.lv+1)) { all.params[i,j] ~ dnorm(0,",1/hypparams[2], ") } } ## All other elements \n", sep = "")) 
		}
	if(num.lv == 1)
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:p) { \n\t\t all.params[i,2] ~ dnorm(0,",1/hypparams[2], ") } \n ", sep = "")) 
			
	if(fc.eff) 
		mod.general.lv <- c(mod.general.lv, paste("\t inter.params ~ dnorm(0,", 1/hypparams[1], ") ## Fourth-corner coefficient \n", sep = ""))
 
	if(site.eff == "fixed") 
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { site.params[i] ~ dnorm(0,", 1/hypparams[1], ") } ## Site fixed effect \n", sep = ""))
	if(site.eff == "ran") 
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { site.params[i] ~ dnorm(0,1/sigma2.site) } \n\t sigma2.site ~ dunif(0,",hypparams[1],") ## Site random effect \n", sep = ""))

	mod.general.lv <- c(mod.general.lv, "\n\t }")
	write(mod.general.lv, file = filename)
	}



## Function to write JAGS code for fitting JSDM fourth corner models
make.jsdmmodel.fourthcorner <- function (family = "binomial", X.eff = "ran", site.eff = "ran", hypparams = c(100,100), fc.eff = FALSE, filename = "jagsjsdmmodel.txt") {
	mod.general.lv <- paste("model { \n\t ## Likelihood \n\t for(i in 1:n) {", sep = "")
	mod.general.lv <- c(mod.general.lv,paste("\t\t Z[i,1:p] ~ dmnorm(eta[i,1:p], Tau[,])", sep = ""))
	mod.general.lv <- c(mod.general.lv,paste("\t\t for(j in 1:p) {", sep = ""))
	
	if(site.eff != "none") site.effect.string <- paste(" + site.params[i]", sep = "")
	if(site.eff == "none") site.effect.string <- NULL

	if(fc.eff) fourthcorner.effect <- paste(" + inter.params*enviro[i]*traits[j]", sep = "")
	if(!fc.eff) fourthcorner.effect <- NULL
		
	if(X.eff == "none") X.effect.string <- NULL
	if(X.eff == "fixed") X.effect.string <- paste(" + enviro[i]*X.params[1] + pow(enviro[i],2)*quad.X.params", sep = "")
	if(X.eff == "ran") X.effect.string <- paste(" + enviro[i]*X.params[j] + pow(enviro[i],2)*quad.X.params", sep = "")

	mod.general.lv <- c(mod.general.lv, paste("\t\t\t eta[i,j] <- all.params[j,1]", site.effect.string, fourthcorner.effect, X.effect.string, sep = ""))

	if(family == "binomial") 			
		mod.general.lv <- c(mod.general.lv, paste("\t\t\t y[i,j] ~ dbern(step(Z[i,j])) } }", sep = ""))
		

	mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:p) { all.params[i,1] ~ dnorm(0,", 1/hypparams[1], ") } ## Species intercept \n", sep = ""))
	if(X.eff != "none") 
		mod.general.lv <- c(mod.general.lv, paste("\t quad.X.params ~ dnorm(0,", 1/hypparams[1], ") \n", sep = ""))
	if(X.eff == "fixed") 
		mod.general.lv <- c(mod.general.lv, paste("\t X.params ~ dnorm(0,", 1/hypparams[1], ") } \n", sep = ""))
	if(X.eff == "ran") {
		mod.general.lv <- c(mod.general.lv, paste("\t sigma2.X.eff ~ dgamma(",1/hypparams[1], ",",1/hypparams[1], ") \n\t mean.X.eff ~ dnorm(0,", 1/hypparams[1], ")", sep = "")) 
		mod.general.lv <- c(mod.general.lv, paste("\t for(j in 1:p) { X.params[j] ~ dnorm(mean.X.eff,1/sigma2.X.eff) } \n", sep = "")) }

	mod.general.lv <- c(mod.general.lv, paste("\t Tau[1:p,1:p] ~ dwish(ID[, ], df) ## Prior for residual covariance matrix \n", sep = ""))
   		
	if(fc.eff) 
		mod.general.lv <- c(mod.general.lv, paste("\t inter.params ~ dnorm(0,", 1/hypparams[1], ") ## Fourth-corner coefficient \n", sep = ""))
 
	if(site.eff == "fixed") 
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { site.params[i] ~ dnorm(0,", 1/hypparams[1], ") } ## Site fixed effect \n", sep = ""))
	if(site.eff == "ran") 
		mod.general.lv <- c(mod.general.lv, paste("\t for(i in 1:n) { site.params[i] ~ dnorm(0,1/sigma2.site) } \n\t sigma2.site ~ dunif(0,",hypparams[1],") ## Site random effect \n", sep = ""))
	
	
	mod.general.lv <- c(mod.general.lv, "\n\t }")
	write(mod.general.lv, file = filename)
	}


## Function for drawing a color gradient legend on plots. Used for unconstrained and residual ordinal plots	
legend.col <- function(col, lev, xadd = 0.1){
 
	opar <- par
	n <- length(col)
	bx <- par("usr")
 
	box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy <- c(bx[3], bx[3])
	box.sy <- (bx[4] - bx[3]) / n
 
	xx <- rep(box.cx, each = 2)
 
	par(xpd = TRUE)
	for(i in 1:n) {
		yy <- c(box.cy[1] + (box.sy * (i - 1)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i - 1)))
		polygon(xx+xadd, yy, col = col[i], border = col[i])
		}

		par(new = TRUE)
		plot(0, 0, type = "n", ylim = c(min(lev), max(lev)), yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
		axis(side = 4, las = 2, tick = FALSE, line = 0, cex.axis = 1)
		par <- opar
	}
	
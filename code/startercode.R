##############
## Starter code for fitting joint models:
## (1) A multivariate GLMM using the lme4 package
## (2) A LVM using the boral package
###############

## We will use the spider dataset available in mvabund, pitfall trap counts of 12 hunting spider species at 28 sites. First load it from mvabund:
library(mvabund)
data(spider)

## Cherry-pick a few variables (lme4 can't handle many responses, especially when we only have 28 sites)
spid.abund = spider$abund[,c(1:3,7,8,12)]

##########
## (1) Fitting a multivariate GLMM using the lme4 package ##
## Use lme4 to fit a multivariate GLMM like in equation 1 of main text, i.e. (random) site effect, and a quadratic effect of soil for each spp. 
## Using a Poisson model with a log link, other options are available via the family argument, as usual.
## Note this function is super-fussy about data - the model didn't converge even for these six species!
##########   

spid.vec = c(as.matrix(spid.abund))
n.site = dim(spid.abund)[1]
n.spp = dim(spid.abund)[2]

## construct a data frame with soil (standardized), species labels, and site labels:
X = data.frame(soil = rep(scale(spider$x[,1]), n.spp), 
               spp = rep(dimnames(spid.abund)[[2]], each=n.site), 
               site = rep(dimnames(spid.abund)[[1]], n.spp) )

## fit the GLMM using lme4 and look at results
library(lme4)
fit.glmm = glmer(spid.vec~0+spp+spp:soil+spp:I(soil^2)+(1|site)+(0+spp|site), data=X, family=poisson())
# Returns some warnings about non-convergence -- the best solution to this is more samples and less variables!!

## The key term in the above is (0+spp|site), which introduces a multivariate random intercept at each site.
## This technique will work for count data but for binary or ordinal data, the variance of the random effect
## needs to be fixed (e.g. at one) for identifiability, which lme4 doesn't do - instead something like MCMCglmm might work.

print(summary(fit.glmm),correlation=FALSE) # To look at estimated parameter values
confint(fit.glmm, method="Wald") # 95% confidence intervals for model parameters. Very approximate though, other more computationally intensive options are available.

## Use corrplot package to calculate and plot residual correlations between species, e.g. possibly due to species interaction etc...
library(corrplot)
vrcorrs=VarCorr(fit.glmm)
corrs=attr(vrcorrs$site.1,"corr")
corrplot(corrs, diag = F, type = "lower", title = "Residual correlations from GLMM", method = "color", tl.srt = 45)



##########
## (2) Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM like in equation 3 of main text, i.e. two latent variables, (random) site effect, and a quadratic effect of soil. 
## Also note that this can fit models to more data (e.g. try full dataset) but takes longer to fit it - for full spdier it will still be a few min, for full aravo dataset over half an hour.
##########   

## Covariates need to be stored as a matrix for boral:
covX <- cbind(scale(spider$x[,1]),scale(spider$x[,1])^2) #soil.dry
## fit the LVM using boral and look at results.  Need version 0.7 or later, available on CRAN.
## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

library(boral)
fit.lvm <- boral(y = spid.abund, 
                 X = covX, 
                 num.lv = 2, 
                 family = "poisson", 
                 row.eff = "random", 
                 save.model = TRUE, 
                 calc.ics = F, 
                 hypparams = c(20,20,20,20))
summary(fit.lvm) # To look at estimated parameter values
fit.lvm$hpdintervals # 95% credible intervals for model parameters.

## compare model coefficients with GLMM:
fixef(fit.glmm)
cbind(fit.lvm$lv.coefs.mean,fit.lvm$X.coefs.mean)
## the slopes are in general agreement, whereas the intercepts differ by quite a bit between the two models.
## The main reason is that in boral, the random intercepts are sampled are a non-zero normal distribution. 

## Dunn-Smyth residual plots  to check model assumption, outliers etc...
plot(fit.lvm)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lvm)

## Use corrplot package to plot residual correlations between species, e.g. possibly due to species interaction etc...
res.cors <- get.residual.cor(fit.lvm)
corrplot(res.cors$cor, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)
## Note correlations are similar to those from the GLMM, the main difference being stronger correlations for the
## two points at the bottom-right of the plot. Adding an extra latent variable to the model would fix this (num.lv=3). 


## Please see help file for the boral function to find further information into the functions available (?boral)
## One nice trick is the get.enviro.cor function, to produce correlations DUE TO environmental variables

### produce my correlation plot

#extract correlation matrices
env.cors <- get.enviro.cor(fit.lvm)
env.tmp<-env.cors$sig.cor
res.cors <- get.residual.cor(fit.lvm)
res.tmp<-res.cors$sig.cor

#make them long
extract_uniquePairDists<-function(dist.mat){

  tmp.df<-data.frame(sp1=NA, sp2=NA, dist=NA)
  i<-0
  j<-0
  for(i in 1:nrow(dist.mat)){
    for(j in 1:ncol(dist.mat)){
      if(upper.tri(dist.mat)[i,j]){
        sp1<-rownames(dist.mat)[i]
        sp2<-colnames(dist.mat)[j]
        dist<-dist.mat[i,j]
        new.row<-data.frame(sp1,sp2,dist)
        tmp.df<-rbind(tmp.df,new.row)
      }

    }
  }

  dist.mat.l<-tmp.df[-1,]
  return(dist.mat.l)

}
View(env.cors)
df.envcors<-extract_uniquePairDists(env.tmp)
df.rescors<-extract_uniquePairDists(res.tmp)

#identify corr abs value, sign
df.cor.list<-list(df.envcors, df.rescors)
data.total.list<-list()
i<-0
for(i in 1:length(df.cor.list)){
  df<-df.cor.list[[i]]
  colnames(df)[3]<-"corr"
  df$absCorr<-abs(df$corr)
  df$sign<-"negative"
  df[df$corr>0,"sign"]<-"positive"
  df[df$corr==0,"sign"]<-"zero"
  data.total.list[[i]]<-df #save this dataframe
}
names(data.total.list)<-c("environment","residual")
library(plyr)
data.total.df<-ldply(data.total.list)

df.summ<-data.total.df
library(tidyr)
df.summ
#df.summ<-gather(df.summ, sign, count, -1)
#df.summ$sign<-mapvalues(df.summ$sign, from=c("numPos","numNeg","numZero"), to=c("positive","negative","zero"))
library(ggplot2)


p<-ggplot(df.summ, aes(x=absCorr, fill=sign)) + 
  geom_histogram(alpha=0.5) + facet_grid(.id~sign) + 
  xlab("Correlation among taxa") + ylab("Frequency")
p



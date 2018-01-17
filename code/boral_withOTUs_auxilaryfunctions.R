
# plot trace

makeMCMCframe <- function(mcmc.obj, mcmc.controls){
  
  # turn the mcmc obj it into a dataframe
  mcmc.df <- data.frame(mcmc.obj[[1]]) # there's only 1 chain for all boral models, hence the 1
  
  # create sample iteration sequence and add it to the df
  sample.iterations <- seq(mcmc.controls$n.burnin+1, mcmc.controls$n.iteration, by = mcmc.controls$n.thin)
  mcmc.df$iteration <- sample.iterations
  
  # gather param cols
  mcmc.df %>%
    gather(key = "param", value = "estimate", -iteration) -> mcmc.df.l
  
  return(mcmc.df.l)

}

plotTrace.lvcoefs <- function(mcmc.obj, mcmc.controls){
  
  mcmc.df.l <- makeMCMCframe(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls)

  mcmc.df.l %>%
    filter(grepl("lv.coefs", param)) %>%
    separate(param, into=c("drop1","drop2","OTUnum","paramNum", "drop3")) %>%
    select(iteration, OTUnum, paramNum, estimate) %>%
    mutate(param = "lv.coefs") -> lv.coefs.df
  
  p <- ggplot(lv.coefs.df, aes(x = iteration, y = estimate, color = OTUnum)) + 
    geom_line() + facet_wrap(~paramNum, scales = "free", nrow = 1) + guides(color = F) +
    ggtitle("LV coefs")
  
  return(p)
  
}

plotTrace.lvs <- function(mcmc.obj, mcmc.controls){
  
  mcmc.df.l <- makeMCMCframe(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls)
  
  mcmc.df.l %>%
    filter(grepl("lvs", param)) %>%
    separate(param, into=c("drop1","OTUnum","paramNum","drop2")) %>%
    select(iteration, OTUnum, paramNum, estimate) %>%
    mutate(param = "lvs") -> lvs.df
  
  p <- ggplot(lvs.df, aes(x = iteration, y = estimate, color = OTUnum)) + 
    geom_line() + facet_wrap(~paramNum, scales = "free") + guides(color = F) +
    ggtitle("LVs")
  
  return(p)
  
}

plotTrace.rowcoefs <- function(mcmc.obj, mcmc.controls){
  
  mcmc.df.l <- makeMCMCframe(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls)
  
  #row coefs 
  #ID 1
  mcmc.df.l %>%
    filter(grepl("row.coefs.ID1", param)) %>%
    separate(param, into=c("drop1","drop2","drop3","paramNum","drop4")) %>%
    select(iteration, paramNum, estimate) %>%
    mutate(param = "rowID1") %>%
    mutate(OTUnum = NA) -> rowid1.df
  #ID2
  mcmc.df.l %>%
    filter(grepl("row.coefs.ID2", param)) %>%
    separate(param, into=c("drop1","drop2","drop3","paramNum","drop4")) %>%
    select(iteration, paramNum, estimate) %>%
    mutate(param = "rowID2") %>%
    mutate(OTUnum = NA) -> rowid2.df
  rowid.df <- full_join(rowid1.df, rowid2.df)
  
  p <- ggplot(rowid.df, aes(x = iteration, y = estimate, color = paramNum)) + 
    geom_line() + facet_wrap(~param, scales = "free") + guides(color = FALSE) +
    ggtitle("rowIDs")
  
  return(p)
  
}

plotTrace.Xcoefs <- function(mcmc.obj, mcmc.controls){
  
  mcmc.df.l <- makeMCMCframe(mcmc.obj = mcmc.obj, mcmc.controls = mcmc.controls)
  
  mcmc.df.l %>%
    filter(grepl("X.coefs", param)) %>%
    separate(param, into=c("drop1", "drop2","OTUnum","paramNum","drop3")) %>%
    select(iteration, OTUnum, paramNum, estimate) %>%
    mutate(param = "Xcoefs") -> Xcoefs.df
  
  p <- ggplot(Xcoefs.df, aes(x = iteration, y = estimate, color = OTUnum)) + 
    geom_line() + facet_wrap(~paramNum, scales = "free") + guides(color = F) +
    ggtitle("Xcoefs")
  
  return(p)
  
}


# shared response and residual correlations

extract_uniquePairDists<-function(dist.mat){
  
  x<-dist.mat
  rowCol <- expand.grid(rownames(x), colnames(x))
  labs <- rowCol[as.vector(upper.tri(x,diag=F)),] #select just the upper triangle without the diagonal no repeats
  df <- cbind(labs, x[upper.tri(x,diag=F)])
  colnames(df) <- c("sp1","sp2","dist")
  
  return(df)
}

extract_cors<-function(corobj, corType){
  
  # Full correlation matrix
  cor <- corobj$cor
  dim(cor)
  cor.df <- extract_uniquePairDists(cor) #make into a long df
  dim(cor.df)
  colnames(cor.df) <- c("otu1","otu2", "cor")
  
  # A correlation matrix containing only the “significant" correlations 
  # whose 95% highest posterior interval does not contain zero. 
  # All non-significant correlations set to zero.
  sig.cor <- corobj$sig.cor
  sig.cor.df <- extract_uniquePairDists(sig.cor) #make into a long df
  colnames(sig.cor.df) <- c("otu1","otu2", "sig.cor")
  
  # join dfs and id 'signif' cells
  cor.df %>%
    left_join(sig.cor.df) %>%
    mutate(signif = ifelse(sig.cor == 0,
                           "no", "yes")) %>%
    select(otu1, otu2, cor, signif) -> df
  
  # replace signif with more specific levels
  df[df$signif == "yes" & df$cor > 0, "signif"] <- "positive"
  df[df$signif == "yes" & df$cor < 0, "signif"] <- "negative"
  
  #update colnames to reflect corType
  colnames(df)[3:4] <- paste(corType, colnames(df)[3:4], sep="_")
  
  return(df)
  
}

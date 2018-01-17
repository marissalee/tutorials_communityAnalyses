Collection of tutorials to analyze community data in R
================

-   [**vegan-centric**](#vegan-centric)
-   [**mvabund-centric**](#mvabund-centric)
-   [**boral-centric**](#boral-centric)

Here are some tutorials that I have come across and have found helpful. I'll also post some code I've been working on. When possible, I've added R scripts and links. Feel free to suggest edits and additions.

**vegan-centric**
-----------------

### Multivariate analyses of ecological communities in R: vegan tutorial

-   Author: Jari Oksanen
-   Type: pdf
-   Source: <http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf>

### Multivariate Statistics

-   Author: Jeff Powell
-   Type: pdf from a Western Sydney HIE course
-   Source: <http://www.hiercourse.com/docs/Rnotes_multivariate.pdf>

### Variation Partitioning

-   Author: Jeff Powell
-   Type: R markdown
-   Source: <http://www.hiercourse.com/docs/variationPartitioning.html>

``` r
#source("code/variationPartitioning.Rmd")
```

### Use distance-based redundancy analysis (vegan::capscale) and model selection (ordistep) to determine which wood traits best explain variation in OTU composition

-   Author: Marissa Lee
-   Type: modified project code
-   Source: <https://github.com/Zanne-Lab/woodEndophytes>

``` r
#source("code/modelselection_capscale.R")
```

**mvabund-centric**
-------------------

### mvabund: an R package for model-based analysis of multivariate abundance data

-   Author: Yi Wang, Ulrike Naumann, Stephen T. Wright, David I. Warton
-   Type: Article in MEE
-   Source: <http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2012.00190.x/full>

### MVABUND package, 2014 updates prepared for Eco-Stats lab

-   Author: David Warton
-   Type: R script from a course
-   Source: <http://eco-stats.blogspot.com.au/2014/09/september-eco-stats-lab-mvabund-package.html>

``` r
#source("code/mvabundTute.R")
```

### Introduction to mvabund

-   Author: Rachel V. Blakey & Andrew Letten
-   Type: Blog post
-   Source: <http://environmentalcomputing.net/introduction-to-mvabund/>

### Use mvabund and model selection to determine which wood traits best explain variation in OTU abundances

-   Author: Marissa Lee
-   Type: modified project code
-   Source: <https://github.com/Zanne-Lab/woodEndophytes>

``` r
#source("code/modelselection_mvabund.R")
```

**boral-centric**
-----------------

### boral – Bayesian Ordination and Regression Analysis of Multivariate Abundance Data in R

-   Author: Francis K.C. Hui
-   Type: Article in MEE
-   Source: <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12514/abstract>
-   Video: <https://www.youtube.com/watch?v=vyMsgyytcUI>

### So Many Variables: Joint Modeling in Community Ecology

-   Author: David I.Warton. F. Guillaume Blanchet, Robert B.O’Hara, Otso Ovaskainen, Sara Taskinen, Steven C.Walker, Francis K.C.Hui
-   Type: Article in TREE
-   Source: <http://www.sciencedirect.com/science/article/pii/S0169534715002402?via%3Dihub>

``` r
#source("code/startercode.R")
#source("code/speciesrichness-analysis.R")
#source("code/fourthcorner-analysis.R")

#which use...
#source("code/fourthcorner-analysis-auxilaryfunctions.R")
#source("code/speciesrichness-auxilaryfunctions.R")
```

### Boral example particularly designed for studies with very diverse communities, e.g. OTU data

-   Author: Marissa Lee
-   Type: modified project code
-   Source: <https://github.com/Zanne-Lab/woodEndophytes>

``` r
#source("code/boral_withOTUs.R")

#which uses...
#source("code/boral_withOTUs_auxilaryfunctions.R")
```

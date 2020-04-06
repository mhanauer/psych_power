---
title: "Diss Power"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
 
Appendix RCode 1: Power for congitive knowledge and skills with normal and non-normal data
```{r}
## Knowledge
pwr.t.test(n = 8, d = 1, type = "paired", alternative = "greater")
pwr.t.test(n = 10, d = .9, type = "paired", alternative = "greater")
pwr.t.test(n = 11, d = .8, type = "paired", alternative = "greater")
pwr.t.test(n = 14, d = .7, type = "paired", alternative = "greater")
pwr.t.test(n = 19, d = .6, type = "paired", alternative = "greater")


### Knowledge not normal
#1 effect size
library(wmwpow)
(.9-.50)/.4
wmwpowd(n = 11, m = 11, distn = "norm(.90,.3)", distm = "norm(.50,.4)", sides = "greater",
alpha = 0.05, nsims=10000)

(.9-.54)/.4
wmwpowd(n = 18, m = 18, distn = "norm(.90,.4)", distm = "norm(.54,.4)", sides = "greater",
alpha = 0.05, nsims=10000)

(.9-.58)/.4
wmwpowd(n = 21, m = 21, distn = "norm(.90,.4)", distm = "norm(.58,.4)", sides = "greater",
alpha = 0.05, nsims=10000)

(.9-.62)/.4
wmwpowd(n = 27, m = 27, distn = "norm(.90,.4)", distm = "norm(.62,.4)", sides = "greater",
alpha = 0.05, nsims=10000)

(.9-.66)/.4
wmwpowd(n = 37, m = 37, distn = "norm(.90,.4)", distm = "norm(.66,.4)", sides = "greater",
alpha = 0.05, nsims=10000)

```


Appendix RCode 2: Replication and power analysis for Cronbach's Alpha

Try to replicate the example in Bonnett (2002)

$$ n_{0} = [8k/(k-1)][Z_{a/2}/ln(e_{1})]^2+2~~ (1)  $$
$$ e_{1} = (1-LL)/(1-UL)~~(2) $$

Replication of Bonnett (2002)
```{r}
k = 4
Za = 1.96
#e1 = (1-.7)/(1-.9)
e1 = 2
n1_first = (8*k)/(k-1)
n1_second = (Za/log(e1))^2
n1 = (n1_first*n1_second)+2
n1

```
Now my study with percision of .7 .9
```{r}
k = 10
Za = 1.96
e1 = (1-.7)/(1-.9)
e1 = 2
n1_first = (8*k)/(k-1)
n1_second = (Za/log(e1))^2
n1 = (n1_first*n1_second)+2
round(n1,0)

round(100/7,2)
```

Appendix RCode 3: Create function based on sim.VSS to create polytomous item structure and demonstrate item loading recovery
```{r}
library(DescTools)
library(psych)
library(StatMeasures)


sim_fifth =  function (ncases = 100, nvariables = 10, nfactors = 1, meanloading = 0.7) 
{
    weight = sqrt(1 - meanloading * meanloading)
    theta = matrix(rnorm(ncases * nfactors), nrow = ncases, ncol = nvariables)
    error = matrix(rnorm(ncases * nvariables), nrow = ncases, 
        ncol = nvariables)
    items = meanloading * theta + weight * error
    items <- apply(items,2, function(x){CutQ(x, breaks = quantile(x, seq(0, 1, by = 0.20)), 
    labels = c(1:5))})
    items = apply(items, 2, function(x){as.numeric(x)})
    return(items)
}


sim_decile =  function (ncases = 1000, nvariables = 10, nfactors = 1, meanloading = 0.5) 
{
    weight = sqrt(1 - meanloading * meanloading)
    theta = matrix(rnorm(ncases * nfactors), nrow = ncases, ncol = nvariables)
    error = matrix(rnorm(ncases * nvariables), nrow = ncases, 
        ncol = nvariables)
    items = meanloading * theta + weight * error
    items <- apply(items, 2,decile)
    return(items)
}

fa

## Turn fifth
dat_fifth = sim_fifth()
# 0.983 at 130
dat_fifth =  sim_fifth(ncases = 120, nvariables = 10)
#dat_fifth
fa_replication  = fa(dat_fifth, 1, rotate="varimax", cor = "poly", correct = 0)
fa_replication
fa_replication$loadings
mean(fa_replication$loadings)



### decile
dat_decile =  sim_decile()
fa_replication  = fa(dat_decile, 1, rotate="varimax", cor = "cor")
fa_replication$loadings
mean(fa_replication$loadings)
decile
```


Ok create simulation
```{r}

n_sample = seq(from = 200, to = 250, by = 10)

efa_power= function(){

n_sample = n_sample
tli_out = list()
rmsea = list()
chi_squre_p = list()
dat_vss = list()
dat_out = list()
for(i in 1:length(n_sample)){
dat_vss[[i]] = sim_fifth(ncases=n_sample[[i]], nvariables=10, nfactors=1, meanloading=.7)
fa_vss[[i]] = fa(dat_vss[[i]], 1, rotate="varimax", cor = "poly", correct = 0)
tli_out[[i]] = fa_vss[[i]]$TLI
rmsea[[i]] = fa_vss[[i]]$RMSEA[1]
chi_squre_p[[i]] = fa_vss[[i]]$PVAL 
dat_out[[i]] = list(tli_out = tli_out[[i]], rmsea = rmsea[[i]], chi_squre_p = chi_squre_p[[i]])
}
return(dat_out)
}
### grab each of them sum them then divide by respective n's
reps = 40
power_efa = replicate(n = reps, efa_power(), simplify = FALSE)
## First 3 tli's are from the first rep meaning they have different sample size.  There are 3 tli's, because there are 3 samples being tested
## the second set of 3 samples is from the second round.  Can we stack them?
power_efa_unlist = round(unlist(power_efa),3)
## Grab the tli for all data sets and reps.  There 
power_efa_matrix = matrix(power_efa_unlist, ncol = 3, byrow = TRUE)
### split data every third row by the number 
colnames(power_efa_matrix) = c("tli", "rmsea", "chi_p")
power_efa_df = data.frame(power_efa_matrix) 
power_efa_df$n = rep(n_sample, reps)
power_efa_df$tli = ifelse(power_efa_df$tli >= .95,1,0)
power_efa_df$rmsea = ifelse(power_efa_df$rmsea <= .05,1,0)
power_efa_df$chi_p = ifelse(power_efa_df$chi_p >= .05,1,0)


power_efa_agg = aggregate(power_efa_df[,1:3], by=list(n=power_efa_df$n), FUN=sum)
# divide by number of reps for %
power_efa_agg[,2:4] =  round(power_efa_agg[,2:4]/reps,3)
power_efa_agg
```
RCode  
```{r}
sim.VSS
dat_vss = sim.VSS(ncases=500, nvariables=10, nfactors=1, meanloading=.7,dichot=TRUE,cut=0)
fa_replication  = fa(dat_vss, 1, rotate="varimax", cor = "tet")
fa_replication$loadings
library(StatMeasures)
scores <- c(1, 4, 7, 10, 15, 21, 25, 27, 32, 35,
            49, 60, 75, 23, 45, 86, 26, 38, 34, 67)

# Create deciles based on the values of the vector
decileScores <- decile(vector = scores)

sim_decile =  function (ncases = 1000, nvariables = 10, nfactors = 1, meanloading = 0.5) 
{
    weight = sqrt(1 - meanloading * meanloading)
    theta = matrix(rnorm(ncases * nfactors), nrow = ncases, ncol = nvariables)
    error = matrix(rnorm(ncases * nvariables), nrow = ncases, 
        ncol = nvariables)
    items = meanloading * theta + weight * error
    items <- apply(items, 2,decile)
    return(items)
}
dat_decile =  sim_decile()

fa_replication  = fa(dat_decile, 1, rotate="varimax", cor = "cor")
fa_replication$loadings
mean(fa_replication$loadings)

```



Try generating a correlation matrix
https://www.rdocumentation.org/packages/clusterGeneration/versions/1.3.4/topics/rcorrmatrix
https://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/
```{r}
# number of observations to simulate
nobs = 100
 
# Using a correlation matrix (let' assume that all variables
# have unit variance 10 items
M = matrix(c(rep(.7,100)), ncol = 10, nrow = 10)
diag(M) = 1
diag(M)
# Cholesky decomposition
L = chol(M)
nvars = dim(L)[1]
nvars 
# R chol function produces an upper triangular version of L
# so we have to transpose it.
# Just to be sure we can have a look at t(L) and the
# product of the Cholesky decomposition by itself
 
t(L)
 
t(L) %*% L
 
 
# Random variables that follow an M correlation matrix
r = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
r = t(r)
rdata = as.data.frame(r)
## Now cut


rdata_cut <- apply(rdata,2, function(x){CutQ(x, breaks = quantile(x, seq(0, 1, by = 0.20)),labels = c(1:5))}) 

### Same for paired or not
cor(mtcars$mpg, mtcars$cyl)
cor.test(mtcars$mpg, mtcars$cyl)
```
Power analysis for pearson correlation
Pearson correlation, two tailed test power of .8 and significance level of .05.
```{r}
library(pwr)
pwr.r.test(n = 100, r = .7)
pwr.r.test

```


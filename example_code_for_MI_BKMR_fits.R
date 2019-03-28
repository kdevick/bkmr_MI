############################################################################
##   Plot code for multiply imputed (MI) datasets using BKMR model fits   ##
##                                                                        ##
##   developed by: Katrina Devick                last updated:  3/28/19   ##
############################################################################


############################################################################
##                           ***** NOTE *****                             ##
## Please read the comments throughout this file and the BKMR_MI_source.R ##
## file. The most important points to highlight are the following:        ##
## 1) All effects throughout this file are calculated as the average      ## 
## change in the outcome for a change in the exposure elements from a     ##
## particular quantile to another quantile calculated across ALL imputed  ##
## datasets. If there are no missing values in the Z matrix (in the       ##
## mixture exposure), this is exactly the same as before.                 ##
## 2) All functions have the option to choose between an "approx" or      ##
## "exact" method. The "exact" method combines the posterior samples      ##
## from all MI fits and uses this posterior chain of #iterations times    ##
## #MI datasets for inference. The "approx" method uses approx estimates  ##
## and std errors from each MI fit and calculates an overall estimate     ##
## and sd using Rubin's 1987 method.                                      ##
## 3) When using the "exact" method, the functions take a while to run,   ##
## so make sure to save the data frames used for plotting and know this   ##
## is to be expected.                                                     ##
############################################################################


## this code assumes you have K BKMR fits and that each of these fits were ran 
## for the same number if MCMC iterations 

## to use this code, first create a list that has each of these
## k=1,...,K fits stored in the kth element 

## fit the K BKMR models
fit1 <- kmbayes(...)
## ...
fitK <- kmbayes(...)


## create a list containing the K BKMR fits
BKMRfits <- list()
BKMRfits[[1]] <- fit1
BKMRfits[[2]] <- fit2
## ...
BKMRfits[[K]] <- fitK



## determine what MCMC iterations of each MI BKMR fit you want to use for inference 
## (you could burnin the first 50% and then keep enough fits to have 1000)
## and then create a vector of the indices to keep

## if the models were fit for 10,000 iterations this would be 
sel.MI <- seq(5001,10000, by = 5)
length(sel.MI) ## double check that this is the correct length (e.g. = 1000)
#### we make the assumption that each of the MI models were fit for the same number of iterations 
#### or that each model was fit for the iterations this vector specifies


## load the source file (change the path to where it is saved)
source("bkmr_MI_source.R")


## load required libraries
library(bkmr)
library(dplyr)
library(magrittr)
library(ggplot2)



########################################################
###     example with multiple simulated datasets     ###
###     (in lieu of MI datasets)                     ###
########################################################



### generate 20 different datasets, fit BMKR and store them in a list 
K <- 20

BKMRfits500 <- list()
for(k in 1:K){
  set.seed(k)
  
  dat <- SimData(n = 500, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X

  fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 10000, verbose = FALSE, varsel = TRUE)
  
  BKMRfits500[[k]] <- fitkm
}
save(BKMRfits500, file="BKMRfits_MI_500.RData")


### load in data for illustration 
load("BKMRfits_MI_500.RData")



######################################################
##         Overall Risk Generation and Plot         ##
######################################################

 
## fit new OverallRiskSummaries function with our MI BKMR fits using approx method
overallrisks.MI <- OverallRiskSummaries.MI(BKMRfits=BKMRfits500, qs = seq(0.1, 0.9, by = 0.05), q.fixed = 0.5, sel = sel.MI, method="approx") 


## rerun the OverallRiskSummaries function but now FIXING the quantile for z2 to 0.25 (for all comparisons)
## NOTE: you can fix one or more elements of the mixture, but the quantile you are fixing them to needs to be the same (q.alwaysfixed)
overallrisks.MI.fixed <- OverallRiskSummaries.MI(BKMRfits=BKMRfits500, qs = seq(0.1, 0.9, by = 0.05), q.fixed = 0.75, q.alwaysfixed = 0.25, index.alwaysfixed = 2, sel = sel.MI, method="approx") 



## When using the exact method, make sure to save this fit since it takes so long to run! 
save(overallrisks.MI, file="overallrisk_MI_500.RData")
## takes approx 9 hours to run for n=500 and 20 MI datasets (fit with simple BKMR model using exact method)
## this will take longer for more complicated BKMR models (model with variable selection, lots
## of covariates, etc.)


## you can now plot this data frame as you normally would for a BKMR fit
ggplot(overallrisks.MI, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +  geom_hline(yintercept=00, linetype="dashed", color="gray")+ 
  geom_pointrange()+ ggtitle("")+ scale_y_continuous(name="estimate")


## plot for the case where you fixed the quantile of 1 (or more) mixture elements 
## if you are making multiple of these plots for different quantiles, you will most likely want to fix 
## the y limits so that the plots are directly comparable 
ggplot(overallrisks.MI.fixed, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) +  geom_hline(yintercept=00, linetype="dashed", color="gray")+ 
  geom_pointrange()+ ggtitle("")+ scale_y_continuous(name="estimate")



######################################################
##     Single Variable Risk Generation and Plot     ##
######################################################


## fit new SingVarRiskSummaries function with our MI BKMR fits
singvarrisk.MI <- SingVarRiskSummaries.MI(BKMRfits=BKMRfits500, qs.diff = c(0.25, 0.75), 
                     q.fixed = c(0.25, 0.50, 0.75), sel=sel.MI, method = "approx")


## dont forget to save this object!
save(singvarrisk.MI, file="singvarrisk_MI_500.RData")
### takes 5.3 hours to run for n=500 and 20 MI datasets (fit with simple BKMR model using exact method)


## rerun the SingVarRiskSummaries.MI function but now FIXING the quantile for z2 to 0.25 (for all comparisons)
## NOTE: you can fix one or more elements of the mixture, but the quantile you are fixing them to needs to be the same (q.alwaysfixed)
## also, you need to specify which.z to NOT include the fixed elements (it will still run if you dont do this, but your plot will look funny)
singvarrisk.MI.fixed <- SingVarRiskSummaries.MI(BKMRfits=BKMRfits500, which.z=c(1,3,4), qs.diff = c(0.25, 0.75), 
                                          q.fixed = c(0.25, 0.50, 0.75), q.alwaysfixed = 0.25, index.alwaysfixed = 2, sel=sel.MI, method = "approx")


## plot the single variable dataframe for the MI fits
ggplot(singvarrisk.MI, aes(variable, est, ymin = est - 1.96*sd,  ymax = est + 1.96*sd, col = q.fixed)) +  geom_hline(aes(yintercept=0), linetype="dashed", color="gray")+ 
  geom_pointrange(position = position_dodge(width = 0.75)) +  coord_flip() + ggtitle("")+ 
  scale_x_discrete(name="Variable")+ scale_y_continuous(name="estimate")



## plot for the case where you fixed the quantile of 1 (or more) mixture elements 
## NOTE: if you do not change which.z to NOT include the fixed element, then the graph will have undesired points
ggplot(singvarrisk.MI.fixed, aes(variable, est, ymin = est - 1.96*sd,  ymax = est + 1.96*sd, col = q.fixed)) +  geom_hline(aes(yintercept=0), linetype="dashed", color="gray")+ 
  geom_pointrange(position = position_dodge(width = 0.75)) +  coord_flip() + ggtitle("")+ 
  scale_x_discrete(name="Variable")+ scale_y_continuous(name="estimate")



######################################################
##       Univariate Risk Generation and Plot        ##
######################################################


## fit new PredictorResponseUnivar function for our MI fits and save the object
## this part can easily be run in parallel by only selecting one z at a time (which.z=1)
## and then rbind-ing the dataframes
univar.MI <- PredictorResponseUnivar.MI(BKMRfits500, ngrid = 50, q.fixed = 0.5, sel = sel.MI, method="approx")
save(univar.MI, file="univar_MI_500.RData")
## takes about 2.7 hours to fun for n=500 and 20 MI datasets (fit with simple BKMR model using exact method)


## plot univariate response functions
ggplot(univar.MI, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + ylab("h(z)") + facet_wrap(~ variable)+ggtitle("")




######################################################
##        Bivariate Risk Generation and Plot        ##
######################################################


## first fit the new PredictorResponseBivar function to our MI fits and save the object 
bivar.MI <- PredictorResponseBivar.MI(BKMRfits = BKMRfits500,  min.plot.dist = 1, sel=seq(5001,10000,by=500), method="approx")
save(bivar.MI, file="bivar_MI_500.RData")
## this take about 9 min to run for n=500 and 20 MI datasets (with simple BKMR model)
## using the approx method 


## Now, apply the UNEDITED function PredictorResponseBivarLevels using the Z matrix
## containing ALL observations from the K MI datasets
Z.MI <- Z.complete.MI(BKMRfits500)
bivar.levels.MI <- PredictorResponseBivarLevels(pred.resp.df = bivar.MI, Z=Z.MI,
                                                          both_pairs = TRUE, qs = c(0.25, 0.5, 0.75))

## create the plot for the bivariate curves
ggplot(bivar.levels.MI, aes(z1, est)) + geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) + ggtitle("h(expos1 | quantiles of expos2)") + xlab("expos1")


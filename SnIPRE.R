#!/usr/bin/Rscript
####
# Run SNIPRE for a given dataset
####

#Example input file:

#        Group.1 PR FR PS FS Total_Syn Trepl nout npop
# NP_001267051.1  3 35  3 41 364.66667  1000   10    8
# NP_001267052.1  9 27  6 12 244.00000  1000   10    8
# XP_003484388.1 10 44  4  1  65.33333  1000   10    8




## Part (1)  Empirical Bayes Implementation  (lme4 package, SnIPRE_source.R)
## Part (2)  Bayesian Implementation (R2WinBUGS package, B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#setwd("/media/data1/bombus/Snipre/snipre_run")
args = commandArgs(trailingOnly = TRUE)
InputFile = args[1]


#################################################################
## Part (1)  Empirical Implementation  (lme4 package)
#################################################################

source("SnIPRE_source.R")
source("my.jags2.R")
library(lme4)
library(R2jags)
library(arm)
library(MASS)

data <- read.table(InputFile, header = TRUE)  # sample data set
data$Trepl =  as.integer(data$Trepl)
data$Tsil =  as.integer(data$Tsil)
snplen = rowSums(data[c(2:5)]) #Number of SNPs counted
len = rowSums(data[c(6,7)]) #length of the Replacement and Silent sites, together
data = (data[(len-snplen)>1,]) #remove any cases where Rep + Syn sites < number of SNPs (odd cases): added by BAH 28-Jul-15 CHECK THESE, they are mistakes!

#SnIPRE <-function(mydata)
# mydata: name of data set;
# mydata must have a header with the following columns: PS, PR, FS, FR, npop, nout, Tsil, Trepl (no particular order)
# outputs 2 objects:  new.dataset & model
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
eb.res = SnIPRE(data)

res = eb.res$new.dataset
model = eb.res$model
write.table(res, file = paste(InputFile, ".ebresults",sep=""), sep  = ",", row.names = FALSE)



#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################

source("B_SnIPRE_source.R")

#data <- read.table(InputFile, header = TRUE)
#BSnIPRE.run <- function(mydata, path = ".", burnin = 500, thin = 5, iter = 2500){
  # path will be where the chains are stored, and must also be where the ".bug" model is located
  # burnin, thin, and iter (number iterations after burnin) are for MCMC samples
BSnIPRE.run(data, burnin = 10000, thin = 4, iter = 15000)

# check to make sure it finished correctly:
# if a "sample" file is in your working directory (getwd()), or the path you sepecified)
# is empty or not there, there is a problem


load("samples")

res.mcmc <- samples

#BSnIPRE <- function(data.mcmc,mydata){
# outputs 2 objects:  new.dataset & effects
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
# the "effects" may be useful if you are interested in estimation
# of population parameters (gamma, constraint) with other assumptions than the PRF

b.res <- BSnIPRE(res.mcmc, data)

bres = b.res$new.dataset

write.table(bres, file = paste(InputFile, ".bayesianresults",sep=""), sep  = ",", row.names = FALSE)


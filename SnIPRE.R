#!/usr/bin/Rscript
####
# Run SNIPRE for a given dataset
####

#Example input file:


args = commandArgs(trailingOnly = TRUE)
InputFile = args[1]


#################################################################
## Mung Data Frame, load datasets
#################################################################
source("B_SnIPRE_source.R")
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


#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################



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
b.res <- BSnIPRE(res.mcmc, data)
bres = b.res$new.dataset
write.table(bres, file = paste(InputFile, ".bayesianresults",sep=""), sep  = ",", row.names = FALSE)


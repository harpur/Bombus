###
#
###













# Load functions and dataframes -------------------------
source("/media/data1/forty3/brock/scripts/VarFunct.r")
gff = read.table(file="/media/data1/bombus/GFF/ref_BIMP_2.0_scaffolds_CDS.txt",header=F)
genelist = read.table(file="/media/data1/bombus/GFF/XP_protein.txt")
gff = gff[gff$V10 %in% genelist$V1,]
bter.snps = read.table(file="/media/data1/bombus/Snipre/snipre_run/FixPoly_Bter.txt",header=F) 
bimp.snps = read.table(file="/media/data1/bombus/Snipre/snipre_run/FixPoly_Bimp.txt",header=F) 
genic.snps = read.table(file="/media/data1/bombus/Snipre/snipre_run/Bter_Bimp_Syn_Poly_Mix.txt",header=F) 
tsil = read.table(file="/media/data1/bombus/Snipre/snipre_run//total_syn_repl_bimp.txt",header=T)
names(tsil)[1]="Group.1"
chrs = gff[c(1,10)]
chrs = chrs[!duplicated(chrs$V10),]; names(chrs)[2]="Group.1"


# Isolate + strand genes and define gene region and utr region -----------------------------
gf.plus = gff[gff$V7=="+",]
gf.plus.st = aggregate(gf.plus$V4, by = list(gf.plus$V10), min)
gf.plus.st$en = aggregate(gf.plus$V5, by = list(gf.plus$V10), max)$x
utr = gf.plus.st$x - 1000 
utr[utr<1] = 1 #if there is less than 1 Kb infront of the gene, set it equal to 1
gf.plus.st$utr = utr

# Isolate - strand genes and define gene region and utr region -----------------------------
gf.neg = gff[gff$V7=="-",]
gf.neg.en = aggregate(gf.neg$V5, by = list(gf.neg$V10), max)
gf.neg.en$en = aggregate(gf.neg$V4, by = list(gf.neg$V10), min)$x
utr = gf.neg.en$x + 1000 
gf.neg.en$utr = utr
gf.utr = rbind(gf.neg.en, gf.plus.st)
#gf.utr$high = apply(gf.utr[c(2,3,4)],1,max) 
#gf.utr$low = apply(gf.utr[c(2,3,4)],1,min)
gf.utr = merge(chrs, gf.utr,by="Group.1")






# Use regions list to get list of Fixed and Polymorphic UTR SNPs -----------------------------
gf.utr$high = apply(gf.utr[c(3,5)],1,max) 
gf.utr$low = apply(gf.utr[c(3,5)],1,min)
bimp.out.snps = c()
for (i in unique(gf.utr$V1)){
	set = gf.utr[gf.utr$V1==i,]
	snps = bimp.snps[bimp.snps$V1==i,]
	blah1=outer(snps$V2,set$high, "<=") 
	blah=outer(snps$V2,set$low, ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T))
	if(is.null(nrow(blah))){
	print(i)
	}else{
	overlap = snps[blah[,1],]
	count = duplicated(overlap$V2) #to ensure I have no overlapping UTRs
	pos = overlap$V2[count=="TRUE"]
	count[overlap$V2 %in% pos] = "TRUE"
	overlap$count = count
	overlap$XP = set[blah[,2],1]
	bimp.out.snps = rbind(bimp.out.snps,overlap)
	}
}



#example of overlapping genes:
#1757816   NT_177094.1  62144  A  T  P  0  TRUE XP_003492399.1
#1757817   NT_177094.1  62157  G  C  P  0  TRUE XP_003492399.1
#1757930   NT_177094.1 117040  A  G  P  0  TRUE XP_012245453.1
#1757816.1 NT_177094.1  62144  A  T  P  0  TRUE XP_012245455.1
#1757817.1 NT_177094.1  62157  G  C  P  0  TRUE XP_012245455.1
#1757930.1 NT_177094.1 117040  A  G  P  0  TRUE XP_012245473.1


bter.out.snps = c()
for (i in unique(gf.utr$V1)){
	set = gf.utr[gf.utr$V1==i,]
	snps = bter.snps[bter.snps$V1==i,]
	blah1=outer(snps$V2,set$high, "<=") 
	blah=outer(snps$V2,set$low, ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T))
	if(is.null(nrow(blah))){
	print(i)
	}else{
	overlap = snps[blah[,1],]
	count = duplicated(overlap$V2) #to ensure I have no overlapping UTRs
	pos = overlap$V2[count=="TRUE"]
	count[overlap$V2 %in% pos] = "TRUE"
	overlap$count = count
	overlap$XP = set[blah[,2],1]
	bter.out.snps =rbind(bter.out.snps,overlap)
	
	}
}


# Get list of overlapping UTRS/Genes -----------------------------
duplic.bimp = as.character(unique(bimp.out.snps$XP[bimp.out.snps$count=="TRUE"])) 
duplic.bter = as.character(unique(bter.out.snps$XP[bter.out.snps$count=="TRUE"]))
duplic = unique(c(duplic.bter, duplic.bimp))


# Output SYN SNPs for every gene and count them -----------------------------
genic.snps = genic.snps[genic.snps$V7=="SYNONYMOUS_CODING",]
genic.snps$V8 = as.character(genic.snps$V8)
snipre.run = (aggregate(genic.snps$V8,by=list(genic.snps$V3),table ))


# Create SNIPRE input file -----------------------------
	#Example:
	#Group PR FR PS FS Tsil Trepl nout npop
	#NP_001267051.1 2 4 3 41 364.666666667 1360.33333333 10 8
	#NP_001267052.1 2 2 6 12 244 872 10 8


#FS and PS
snipre.run = (aggregate(genic.snps$V8,by=list(genic.snps$V3),function(x) length(x[x=="FS"]) ))
names(snipre.run)[2]="FS"
snipre.run$PS = (aggregate(genic.snps$V8,by=list(genic.snps$V3),function(x) length(x[x=="PS"]) ))$x

#FR and PR
bter = (aggregate(bter.out.snps$V5,by=list(bter.out.snps$XP),function(x) length(x[x=="F"]) ))
bter$PR = (aggregate(bter.out.snps$V5,by=list(bter.out.snps$XP),function(x) length(x[x=="P"]) ))$x
names(bter)[2] = "FR"
bimp = (aggregate(bimp.out.snps$V5,by=list(bimp.out.snps$XP),function(x) length(x[x=="P"]) ))
names(bimp)[2] = "PR"

#assemble the file
snipre.run = merge(snipre.run, bter, by = "Group.1",all.x=T)
snipre.run = merge(snipre.run, bimp, by = "Group.1",all.x=T)
snipre.run$PR = snipre.run$PR.x + snipre.run$PR.y; snipre.run$PR.x = snipre.run$PR.y = NULL
snipre.run = merge(snipre.run, tsil, by = "Group.1",all.x=T)	
snipre.run = merge(snipre.run, gf.utr, by = "Group.1",all.x=T)
snipre.run$Trepl = as.numeric(snipre.run$high) - as.numeric(snipre.run$low)+1
snipre.run = snipre.run[c(1,5,4,3,2,6,15)]
snipre.run$nout = rep(10, nrow(snipre.run))
snipre.run$npop = snipre.run$nout-2
snipre.run$PR[is.na(snipre.run$PR)]	= 0
snipre.run$FR[is.na(snipre.run$FR)]	= 0
snipre.run = snipre.run[snipre.run$Trepl>1,]
names(snipre.run)[6] = "Tsil"
#remove overlapping genes/UTRs
snipre.run = snipre.run[!(snipre.run$Group.1 %in% duplic),] # 7519 included because they don't overlap
write.list(snipre.run, file="/media/data1/bombus/Snipre/snipre_run/UTRSNIPREnooverlapping.input")



# Running SNIPRE -----------------------------
	#See SnIPRE_example_script.R

#Rscript SnIPRE.R /media/data1/bombus/Snipre/snipre_run/UTRSNIPRE.input
	#write.table(bres, file = paste("UTR", ".bayesianresults",sep=""), sep  = ",", row.names = FALSE)


	
	
	
#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################

source("/media/data1/bombus/Snipre/snipre_doc/B_SnIPRE_source.R")
source("/media/data1/bombus/Snipre/snipre_doc/my.jags2.R")
library(lme4)
library(R2jags)
library(arm)


data <- snipre.run#read.table(InputFile, header = TRUE)
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

b.res.overlap <- BSnIPRE(res.mcmc, data)

bres.overlap = b.res.overlap$new.dataset
	

	hist(bres.overlap$BSnIPRE.gamma,breaks=50)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


	

#As of now, I am NOT altering for coverage and I AM altering for overlapping SNPs.






# Get out Trep	 -----------------------------
	#Trep is the number of sites (/1000) covered with at least 2 reads in both species
	#I need that number AND the sites so I can filter out SNPs.



#genelist filter.....






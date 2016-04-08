#!/usr/bin/Rscript
###
# Count Fixed and Polymorphic Sites in Haploid Data
###









#Creates a SNIPRE input table using SNPEFF outputs and merged VCF file as inputs as well as sample counts
	#Example (not real data) here:
# GeneID PR FR PS FS Tsil Trepl nout npop
# NP_001267051.1  3 35  3 41 364.66667  1000   10    8
# NP_001267052.1  9 27  6 12 244.00000  1000   10    8
# XP_003484388.1 10 44  4  1  65.33333  1000   10    8

#Arguments ----------------------------------------------
i = "bombus.vcf" #argument here
n1 = 10 #args, sample size for pop 1, should be first samples in list
n2 = 8 #args, sample size for pop 2
snpeff  = "exons.eff" #argument here


# Load functions ------------------------------------------
source("/git/GenomeR/VCFFunctions.r")

mk.test<-function(x){
	mat = matrix(nr=2,nc=2)
	mat[,1] = c(x[4],x[3])
	mat[,2] = c(x[2],x[1])
	fish = fisher.test(mat)$p.value
	ni = (x[1] / x[3]) / (x[2] / x[4])
	return(c(fish,ni))
}

# Load VCF and Putative Functional Roles -------------------------------------------------
vcf = Read.VCF()
snpeff  = read.table(file=snpeff, header=F) 
snpeff  = snpeff[c(1,2,9,13)];snpeff$ID = paste(snpeff$V1, snpeff$V2, sep="_")
sites  = read.table(file="/data2/bombus/git/data/total_syn_repl_bimp.txt", header=T) 

# Identify sample sizes and identities ----------------------
n1.last = 10 + n1 -1
n2.first = n1.last + 1
n2.last = n2.first + n2 -1

# Trim trialleles ------------------------------------------
alleles = sapply(vcf,function(x) nchar(x[5]))
vcf = vcf[-which(alleles>1)]

#Get list of all sites --------------------------------------
chrom = sapply(vcf,function(x) return(c(x[1],x[2]))) #list containing  alls caffolds and positions

#Count number of polymorphic sites and missing sites for N1 ------------------
ones = sapply(vcf, function(x) (length(x[10:n1.last][grep("^1:", x[10:n1.last])]))) #counts number of "1:'s"
zers = sapply(vcf, function(x) (length(x[10:n1.last][grep("^0:", x[10:n1.last])]))) #counts number of "0:'s"
fixed = sapply(vcf, function(x) (length(x[10:n1.last][grep("^[.]", x[10:n1.last])]))) #counts number of ".'s"
miss = sapply(vcf, function(x) (length(x[10:n1.last][grep("^[.]:", x[10:n1.last])]))) #counts missing data

#Assess whether Fixed or Poly in  N1 ----------------------
size = length(10:n1.last)
size = rep(size, length(miss))
size = size - miss

FP.n1 = ones
FP.n1[which(fixed==size)] ="NA"
FP.n1[which(ones==size)] ="F"
FP.n1[which(ones!=size & (zers+ones) == size)] ="P"


#Count number of polymorphic sites and missing sites for N2 ------------------
ones.n2 = sapply(vcf, function(x) (length(x[n2.first:n2.last][grep("^1:", x[n2.first:n2.last])]))) #counts number of "1:'s"
zers.n2 = sapply(vcf, function(x) (length(x[n2.first:n2.last][grep("^0:", x[n2.first:n2.last])]))) #counts number of "0:'s"
fixed.n2 = sapply(vcf, function(x) (length(x[n2.first:n2.last][grep("^[.]", x[n2.first:n2.last])]))) #counts number of ".'s"
miss.n2 = sapply(vcf, function(x) (length(x[n2.first:n2.last][grep("^[.]:", x[n2.first:n2.last])]))) #counts missing data

#Assess whether Fixed or Poly in  N2 ----------------------
size.n2 = length(n2.first:n2.last)
size.n2 = rep(size.n2, length(miss.n2))
size.n2 = size.n2 - miss.n2

FP.n2 = ones.n2
FP.n2[which(fixed.n2==size.n2)] ="NA"
FP.n2[which(ones.n2==size.n2)] ="F"
FP.n2[which(ones.n2!=size.n2 & (zers.n2+ones.n2) == size.n2)] ="P"


#Compare both and identify fixations and Polymorphisms across the lineage ----------------------
FP = rep(0, length(FP.n1))
FP[which(FP.n2=="NA")] = FP.n1[which(FP.n2=="NA")]
FP[which(FP.n1=="NA")] = FP.n2[which(FP.n1=="NA")]
FP[which(FP.n1!="NA" & FP.n2!="NA" & FP.n1!=FP.n2)] = "P" #these are "FP" cases
FP[which(FP.n1!="NA" & FP.n2!="NA" & FP.n1==FP.n2)] =  FP.n1[which(FP.n1!="NA" & FP.n2!="NA" & FP.n1==FP.n2)] 

#Create fixed and polymorphic data frame  ----------------------
FP = data.frame(cbind(chrom = chrom[1,], pos = chrom[2,], FP))
FP$ID = paste(FP$chrom, FP$pos, sep="_")

#Mung SNPEFF dataframe  ----------------------
snpeff = snpeff[which(snpeff$ID %in% FP$ID),] #remove tri-alleles
n_occur = data.frame(table(snpeff$ID))
overlapping.genes = unique(as.character(snpeff$V9[which(snpeff$ID %in% n_occur$Var1[n_occur$Freq > 1])]))
write.table(overlapping.genes, file="OverlappingGenes", col.names=F, row.names=F, quote=F)
snpeff = snpeff[-which(snpeff$ID %in% n_occur$Var1[n_occur$Freq > 1]),] # Remove overlapping genes (this removes 16 genes , 48 SNPs) 

#Update SYN calls   ----------------------
	#"NON_SYNONYMOUS_START" and "SYNONYMOUS_STOP", by SNPEFF definition are both synonymous
nsyn.update = as.character(snpeff$V13)
nsyn.update[nsyn.update=="NON_SYNONYMOUS_START"] ="SYNONYMOUS_CODING"
nsyn.update[nsyn.update=="SYNONYMOUS_STOP"] ="SYNONYMOUS_CODING"
snpeff$V13 = nsyn.update

#Merge data frames, then build MK table ----------------------
MK =  merge(snpeff, FP, by = "ID")
MK$MK = paste(MK$V13, MK$FP, sep="_")
MK  = aggregate(MK$MK, by = list(MK$V9), table)
write.table(MK, file="MKtableraw")
MK = read.table(file="MKtableraw",header=T)
names(MK) = c("GeneID","FR","PR","FS", "PS")
MK = MK[c(1,3,2,5,4)]


# Build Snipre dataframe -----------------------------------
MK.snipre = merge(MK, sites, by="GeneID")
MK.snipre = MK.snipre[-c(8)]
MK.snipre$nout = rep(n2,nrow(MK.snipre))
MK.snipre$npop = rep(n1,nrow(MK.snipre))
write.table(MK.snipre, file="MK.snipre", col.names = T, row.names=F, quote = F)

# Estimate alpha, NI, MK test -----------------------------------
MK.snipre1 = MK.snipre[c(2,3,4,5)] 
mk.res = apply(MK.snipre1, 1, function(x) mk.test(x))

fish = mk.res[1,]
ni = mk.res[2,]
ni[is.na(ni)] = 0
ni[ni=="-Inf"] = 0
ni[ni=="Inf"] = 0
MK.snipre$NI = ni 
MK.snipre$MKp = fish

write.table(MK.snipre, file="MKresults", col.names = T, row.names=F, quote = F)















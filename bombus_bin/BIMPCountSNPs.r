###
#
###


#after SNPEff_extractSNPs.sh
#Counts number of unique SNPs in exons, introns, intergenic




# Load Depth Information --------------------------------
depth = read.table(file="out.ldepth",header=T)
depth$SNP = with(depth, paste(CHROM, POS, sep="_"))
depth$CHROM = NULL
depth$POS = NULL

# Load Exons --------------------------------
exons = read.table(file="exons.filt.1.eff",header=F)
exons$SNP = with(exons, paste(V1, V2, sep="_"))
exons = merge(exons, depth, by = "SNP")
exons$type = rep("exons", nrow(exons))

# Load  Upstream (5kB) --------------------------------
upstr = read.table(file="upstr.eff",header=F)
upstr$SNP = with(upstr, paste(V1, V2, sep="_"))
	#upstr$V12 contains the length away fromt he start of a gene a SNP is. I'll consider only 5Kb
upstr = merge(upstr, depth, by = "SNP")
upstr$type = rep("upstr", nrow(upstr ))	
	
	
# Load Intergenic --------------------------------
interg = read.table(file="interg.eff",header=F)
interg$SNP = with(interg, paste(V1, V2, sep="_"))
interg=interg[!(interg$SNP %in% upstr$SNP),]
interg = merge(interg, depth, by = "SNP")
interg$type = rep("interg", nrow(interg))


# Load Introns --------------------------------
introns = read.table(file="intron.eff",header=F)
introns$SNP = with(introns, paste(V1, V2, sep="_"))
introns=introns[!(introns$SNP %in% upstr$SNP),]
introns = merge(introns, depth, by = "SNP")
introns$type = rep("introns", nrow(introns))


# Put everything together --------------------------------
snps = rbind(introns[c(1,13,15)], exons[c(1,19,21)], interg[c(1, 10,12)])
boxplot(snps$SUM_DEPTH~snps$type)











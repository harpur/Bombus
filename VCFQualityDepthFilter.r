#!/usr/bin/Rscript
###
# Trim VF by quality and depth
###




args <- commandArgs(trailingOnly = TRUE)
fil = as.character(args[1]) #file name 

#Output site depth and quality with vcftools --------
system(paste("vcftools --vcf ",fil,"--site-mean-depth", sep=" "))
system(paste("vcftools --vcf ",fil,"--site-quality", sep=" "))

#get max and min depth values ---------------------
depth=read.table(file="out.ldepth.mean",header=T)
maxdp = 1.5*IQR(depth$MEAN_DEPTH)+quantile(depth$MEAN_DEPTH,0.75) #max depth outlier 
mindp = quantile(depth$MEAN_DEPTH,0.25)-1.5*IQR(depth$MEAN_DEPTH) #min depth outlier

qual=read.table(file="out.lqual",header=T)
maxq = 1.5*IQR(qual$QUAL)+quantile(qual$QUAL,0.75) #max qual outlier 
minq = quantile(qual$QUAL,0.05) #min qual outlier
qual = qual[which(qual$QUAL< minq | qual$QUAL > maxq),]	

#Trim low and high depth sites with vcftools --------
system(paste("vcftools --vcf", fil," --min-meanDP", mindp, "--max-meanDP", maxdp, "--recode --out out.indel.dp", sep=" "))

#Remove Quality outliers ----------------------------
write.table(qual[c(1,2)], file="QualityOutlierSNPs",col.names=F,row.names=F,quote=F)
system("vcftools --vcf out.indel.dp.recode.vcf --exclude-positions QualityOutlierSNPs --recode --remove-filtered-all  --out out.indel.dp.q")
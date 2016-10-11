###
# File Prep for Bombus Balancing Scans
###


#BTER:
#1) Use GATK UnifiedGenotype - this takes ~19hrs hrs using 5+15cpus.

#SNPs
gatk -R /media/data1/bombus/BterrFullGenome.fa -T UnifiedGenotyper -I Bter-87_001.sorted.gatk.bam  -I Bter-90_001.sorted.gatk.bam -I Bter-91_001.sorted.gatk.bam -I Bter-81_001.sorted.gatk.bam -I Bter-82_001.sorted.gatk.bam -I Bter-83_001.sorted.gatk.bam -I Bter-84_001.sorted.gatk.bam -I Bter-85_001.sorted.gatk.bam -o Bter.raw.vcf -stand_call_conf 60.0 -stand_emit_conf 40.0 -dcov 200 --min_base_quality_score 20 -nt 8 -glm SNP -ploidy 1
  #removed 88 and 89, they are too low quality and highly hetero
	#-I Bter-88_001.sorted.gatk.bam -I Bter-89_001.sorted.gatk.bam
  
  
  
#INDELS
gatk -R /media/data1/bombus/BterrFullGenome.fa -T UnifiedGenotyper -I Bter-87_001.sorted.gatk.bam -I Bter-88_001.sorted.gatk.bam -I Bter-89_001.sorted.gatk.bam -I Bter-90_001.sorted.gatk.bam -I Bter-91_001.sorted.gatk.bam -I Bter-81_001.sorted.gatk.bam -I Bter-82_001.sorted.gatk.bam -I Bter-83_001.sorted.gatk.bam -I Bter-84_001.sorted.gatk.bam -I Bter-85_001.sorted.gatk.bam -o Bter.raw.indel.vcf -stand_call_conf 60.0 -stand_emit_conf 40.0 -dcov 200 --min_base_quality_score 20 -nt 12 -glm INDEL



##
#Make first filtered file (RAW) called Bter.raw.recode.vcf
#gatk -R /media/data1/bombus/BterrFullGenome.fa -T VariantFiltration -V  Bter.raw1.vcf --filterExpression "QD < 5.0 || FS > 40.0 || MQ < 25.0 || DP < 100.0"   --filterName "mpileFilters" -o  Bter.raw2.vcf
#
vcftools --vcf Bter.raw2.vcf --recode --remove-filtered-all --out Bter.raw



##
#Make second Filtered File (INDELS) Bter.indel.recode.vcf
gatk -R /media/data1/bombus/BterrFullGenome.fa -T VariantFiltration -V  Bter.raw.recode.vcf -o  Bter.indel.vcf --mask Bter.raw.indel.vcf --maskExtension 10 --maskName "InDel" 

vcftools --vcf Bter.indel.vcf --recode --remove-filtered-all --out Bter.indel


##
#Make third file (depth), called Bter.indel.dp.recode.vcf
vcftools --vcf Bimp.indel.recode.vcf --site-mean-depth
#Coded in a separate Rscript, but here is some of the code within:
R
depth=read.table(file="out.ldepth.mean",header=T)
maxdp = 1.5*IQR(depth$MEAN_DEPTH)+quantile(depth$MEAN_DEPTH,0.75) #maxdp
mindp = quantile(depth$MEAN_DEPTH,0.25)-1.5*IQR(depth$MEAN_DEPTH) #mindp
d1 = depth[depth$MEAN_DEPTH <= mindp,]
d2 = rbind(d1, depth[depth$MEAN_DEPTH >= maxdp,])
write.table(d2[c(1,2)], file="BIMPDepthSNPs",col.names=F,row.names=F,quote=F)

vcftools --vcf Bimp.indel.dp.recode.vcf  --exclude-positions BIMPDepthSNPs --recode --remove-filtered-all  --out Bter.indel.dp


##
#Make forth file (BLAST filter) Bter.indel.dp.recode.vcf
vcftools --vcf Bimp.indel.dp.recode.vcf  --exclude-positions /media/data1/bombus/MaskedSNPs/BimpAllSNPs --recode --remove-filtered-all  --out Bimp.indel.dp.bl



##
#Make forth file (diploid drone calls) Bter.indel.dp.bl.dd.recode.vcf
	#Bter.raw.DIPLOID.vcf contains a set of SNOs called when diploid
	#Had to remove 88 and 89
	#vcftools --vcf Bter.raw.DIPLOID.vcf  --remove-indv "Bter-89" --remove-indv "Bter-88" --recode

vcftools --vcf Bter.raw.DIPLOID.vcf  --hardy 

R
hwe=read.table(file="out.hwe",header=T)
hets=apply(hwe[3],1,function(x) unlist(strsplit(x, "/"))[2])
hwe1=hwe[c(1,2)][which(hets>0),]
write.table(file="/media/data1/bombus/MaskedSNPs/DiploidBterr",hwe1,col.names=F,row.names=F,quote=F)


vcftools --vcf Bimp.indel.dp.bl.recode.vcf   --exclude-positions /media/data1/bombus/MaskedSNPs/DiploidBimp --recode --remove-filtered-all --out Bimp.indel.dp.bl.dd




#There it is, I get five files:
	#Bter.raw.recode.vcf
	#Bter.indel.recode.vcf
	#Bter.indel.dp.recode.vcf
	#Bter.indel.dp.recode.vcf
	#Bter.indel.dp.bl.dd.recode.vcf
















#Calling NSYN and SYN SNPs:
#http://snpeff.sourceforge.net/SnpEff_manual.html
	#search "databases"
	#Update the .config file 

java -jar /usr/local/lib/snpEff2/snpEff.jar build -gff3 -v Bterr


#Run SNPEff
	#Run on AfrSNPs.raw.vcf, it contains all SNPs 
java -jar /usr/local/lib/snpEff2/snpEff.jar Bterr -o txt Bter.raw.recode.vcf  -no-downstream -no-upstream  > Bterr.snpeff.eff

	#Will NEED to remove genes with premature stop codons and any trimmed genes in Daria's analyses

sed '/SYNONYMOUS/ !d' Bterr.snpeff.eff > exons.eff
sed '/WARNING/ !d' exons.eff > warnexons.eff
sed '/WARNING/ d' exons.eff > exons1.eff #renamed exons.eff



#Make a vcf file of ONLY the NSYN SNPs and ONLY the SYN SNPs:
	#wrote out NSYN.snp and SYN.snp from R using exons1.eff
	nsyn=read.table(file="/media/data1/forty3/brock/align/exons.eff")
	nsyn=nsyn[,-c(5:7,9,10,11,12,16,17)]
	nsyn$nsyn=gsub("_.*","",nsyn$V13)
	nsyn$V13=NULL
	syn=nsyn[nsyn$nsyn!="NON",]
	write.table(syn[c(1,2)],file="SYN.snp",quote=F,col.names=F,row.names=F)
	nsyn=nsyn[nsyn$nsyn=="NON",]
	write.table(nsyn[c(1,2)],file="NSYN.snp",quote=F,col.names=F,row.names=F)

	nsyn=read.table(file="/media/data1/forty3/brock/align/warnexons.eff")
	nsyn=nsyn[,-c(5:8,10,11,12,13,16,17,18)]
	nsyn$nsyn=gsub("_.*","",nsyn$V14)
	syn=nsyn[nsyn$nsyn!="NON",]
	write.table(syn[c(1,2)],file="SYN.snp",quote=F,col.names=F,row.names=F,append=T)
	nsyn=nsyn[nsyn$nsyn=="NON",]
	write.table(nsyn[c(1,2)],file="NSYN.snp",quote=F,col.names=F,row.names=F,append=T)




 
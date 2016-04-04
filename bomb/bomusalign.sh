###
# Call SNPs for Bombus Genomes 
###



#Environment variables -------------------------------
export BOMBUS=/media/data2/bombus/
export MASK=/media/data1/bombus/MaskedSNPs/



			export BAM=/media/data1/afz/BAM

			export VCF=/media/data1/afz/VCF
			export AFZ=/media/data1/afz
			export GIT=/media/data1/git
			export FASTQ=/media/data1/fastq
			export REF=/media/data1/fastq


			/media/data1/bombus/BimpFullGenome.fa



# Call SNPs with GATK -----------------------------------
	#All 
ls *.bam > bams.list  #assemble list of bams in species directory

#Call SNPs
gatk -R $BOMBUS/BimpFullGenome.fa -T UnifiedGenotyper \
	-I bams.list  \
	-o out.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 20 -glm SNP  \
	-ploidy 1 &
 
#Call INDELs
gatk -R $BOMBUS/BimpFullGenome.fa -T UnifiedGenotyper \
	-I bams.list  \
	-o out.raw.indels.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 20 -glm indel  \
	-ploidy 1 &
 

#Call SNPs - assume diploid
gatk -R $BOMBUS/BimpFullGenome.fa -T UnifiedGenotyper \
	-I bams.list  \
	-o out.het.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 20 -glm SNP  \
	-ploidy 2
		
		SCRIPT:
		R
			hwe=read.table(file="out.hwe",header=T)
			hets=apply(hwe[3],1,function(x) unlist(strsplit(x, "/"))[2])
			hwe1=hwe[c(1,2)][which(hets>0),]
			write.table(file="/media/data1/bombus/MaskedSNPs/DiploidBterr",hwe1,col.names=F,row.names=F,quote=F)


			vcftools --vcf Bimp.indel.dp.bl.recode.vcf   --exclude-positions /media/data1/bombus/MaskedSNPs/DiploidBimp --recode --remove-filtered-all --out Bimp.indel.dp.bl.dd
			 
 
 # Filter Variants  -----------------------------------
gatk -R $BOMBUS/BimpFullGenome.fa -T VariantFiltration 
	-V  out.raw.vcf \
	--filterExpression "QD < 5.0 || FS > 40.0 || MQ < 25.0 "   \
	--filterName "mpileFilters" \
	--mask out.raw.indels.vcf  \
	--maskExtension 10 \
	--maskName "InDel" \
	-o  out.vcf

vcftools --vcf out.vcf --recode --remove-filtered-all --out out.indel #will be called out.indel.recode.vcf
 
  

# Filter Depth outliers  -----------------------------------
vcftools --vcfout.indel.recode.vcf --site-mean-depth

Rscript  
	depth=read.table(file="out.ldepth.mean",header=T)
	maxdp = 1.5*IQR(depth$MEAN_DEPTH)+quantile(depth$MEAN_DEPTH,0.75) #maxdp
	mindp = quantile(depth$MEAN_DEPTH,0.25)-1.5*IQR(depth$MEAN_DEPTH) #mindp
	d1 = depth[depth$MEAN_DEPTH <= mindp,]
	d2 = rbind(d1, depth[depth$MEAN_DEPTH >= maxdp,])
	write.table(d2[c(1,2)], file="BIMPDepthSNPs",col.names=F,row.names=F,quote=F)

vcftools --vcf out.recode.vcf  --exclude-positions BIMPDepthSNPs --recode --remove-filtered-all  --out out.indel.dp.

 
# Filter by BLAST ----------------------------------- MAYBE!
 vcftools --vcf Bimp.indel.dp.recode.vcf  --exclude-positions $MASK/BimpAllSNPs --recode --remove-filtered-all  --out Bimp.indel.dp.bl
 
 

# output final VCF file ----------------------------------- 








#Calling synonymous and non-synonymous sites ------------------



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













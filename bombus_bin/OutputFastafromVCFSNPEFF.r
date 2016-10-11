###
# Extract gene regions of a target gene and output genotypes for all individuals in a VCF file
###

#Useage
#requirements: Haploid Individual










#XP_01224403.1 #no SNPs




# Functions --------------------------------
source("/media/data1/forty3/brock/scripts/ReferencetoDataframe.r")
source("/media/data1/forty3/brock/scripts/VCFFunctions.r")



# Get out SNPs in and near genes of interest -------------------------------
	#Automate this later....
system("sed '/XP_012247565.1/ !d' /media/data1/bombus/snpeff_Bimp/Bimp.filt.1.snpeff.eff > targ.gene.eff") #output SNPs around target gene using SNPEFF
system("cut -f -2 targ.gene.eff > output_filename") #Output a VCF trim file
system("vcftools --vcf /media/data1/bombus/vcf_Bimp/Bimp.indel.dp.bl.recode.vcf --positions output_filename --recode --out targgene") #create trimmed VCF file


#BTER!
# Get out SNPs in and near genes of interest -------------------------------
	#Automate this later....
system("sed '/XP_012247565.1/ !d' /media/data1/bombus/snpeff_Bimp/Bimp.filt.1.snpeff.eff > targ.gene.eff") #output SNPs around target gene using SNPEFF
system("cut -f -2 targ.gene.eff > output_filename") #Output a VCF trim file
system("vcftools --vcf /media/data1/bombus/vcf_Bimp/Bimp.indel.dp.bl.recode.vcf --positions output_filename --recode --out targgene") #create trimmed VCF file







# Load Genome ----------------------------------
AMEL = Reftodf(data.frame = "/media/data1/bombus/BimpFullGenome.fa" , nms = "N")

# Load VCF File  ----------------------------------
i = "targgene.recode.vcf"
vcf = Read.VCF()

# Extract Bee IDs ------------------
ids = unlist(vcf[1])[c(10:length(unlist(vcf[1])))] 

# Extract and output sequences in .fasta --------------
for(i in 1:length(ids)){

	# Extract SNPs and their positions ------------------
	pos = unlist(lapply(vcf,function(x)  unlist(x[2])))
	pos = as.numeric(pos[-1])
	alt.raw = as.character(unlist(lapply(vcf,function(x)  unlist(x[5]))))
	alt.raw = alt.raw[-1]
	alt = substr(alt.raw,1,1)
	alt2 = substr(alt.raw,3,3) #for tri-alleles
	ref = as.character(unlist(lapply(vcf,function(x)  unlist(x[4]))))
	ref = ref[-1]

	# extract bee genotype ---------------------
	geno = as.character(unlist(lapply(vcf,function(x)  unlist(x[i + 9])))) 
	geno = geno[-1]
	geno = gsub(":.*","",geno)
	geno[which(geno=="0")] = ref[which(geno=="0")]
	geno[which(geno=="1")] = alt[which(geno=="1")]
	geno[which(geno==".")] = "N"
	geno[which(geno=="2")] = alt2[which(geno=="2")]
	
	# Define and Extract Gene Region ------------------
	chr = unlist(vcf[2])[1]
	st = as.numeric(unlist(vcf[2])[2])
	en = as.numeric(unlist(vcf[length(vcf)])[2])
	sequ = as.character(AMEL$seqs[AMEL$chrom==chr])
	sequ = unlist(strsplit(sequ, split=""))

	# Input genotypes into sequence ---------------
	sequ[pos] = geno
	sequ = sequ[c(st:en)]
	sequ = paste(sequ, collapse="")
	
	# Write out .fasta ------------------------
	sink(file=paste("GENEID", ".fas",sep=""), append=T)
	cat(c(">",ids[i]));cat("\n")
	cat(sequ);cat("\n")		
	sink()
		
}
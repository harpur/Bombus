###
# Call SNPs for Bombus Genomes 
###



#Environment variables -------------------------------
export BOMBUS=/media/data2/bombus
export MASK=/media/data1/bombus/MaskedSNPs
export SNPEFF=/usr/local/lib/snpEff2


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
	
#Output hetero SNPs --------------------
vcftools --vcf out.het.raw.vcf --hardy #2792802 out.hwe
sed '/[/]0[/]/ d' out.hwe > test.out
sed '/[/]1[/]/ d' test.out > test1.out #447528SNPs 
cut -f 1,2 test1.out > putativeCNV.list 
rm test.out
rm test1.out
grep -Fwf putativeCNV.list out.raw.vcf > cnv.vcf
sed '/^#/ !d' out.raw.vcf > VCFheader
cat VCFheader cnv.vcf > cnv1.vcf
rm cnv.vcf
mv cnv1.vcf cnv.vcf

# Filter Variants  -----------------------------------
gatk -R $BOMBUS/BimpFullGenome.fa -T VariantFiltration \
	-V out.raw.vcf \
	--filterExpression "QD < 5.0 || FS > 40.0 || MQ < 25.0 "   \
	--filterName "mpileFilters" \
	--mask out.raw.indels.vcf  \
	--maskExtension 10 \
	--maskName "InDel" -o  out.vcf 

gatk -R $BOMBUS/BimpFullGenome.fa -T VariantFiltration \
	-V out.vcf \
	--mask cnv.vcf \
	--maskExtension 5 \
	--maskName "cnv" \
	-o  out2.vcf
rm out.vcf 
mv out2.vcf out.vcf #2488716

vcftools --vcf out.vcf --recode --remove-filtered-all --max-alleles 2 --out out.indel #will be called out.indel.recode.vcf
#wc -l 	out.indel.recode.vcf 2192172 (88% of original SNPs in BIMP)

# Filter Depth and Quality outliers  -----------------------------------
VCFQualityDepthFilter.r "out.indel.recode.vcf" 

# Calling synonymous and non-synonymous sites -------------------------
	#note, v4.0 of SNPEFF no longer supports TXT format (I think)...annoying.
java -jar $SNPEFF/snpEff.jar Bimp \
	-o txt \
	out.indel.dp.q.recode.vcf \
	-no-downstream \
	-no-upstream > out.snpeff.eff
	
sed '/SYNONYMOUS/ !d' out.snpeff.eff > exons.eff #take only SNPS within exons from output file
sed '/WARNING/ !d' exons.eff > warnexons.eff #take out warnings from output file
sed -i '/WARNING/ d' exons.eff 


# make VCF file of only SYN and NSYN SNPs -------------------------
cut -f 1,2 exons.eff  > synnsyn.list 
vcftools --vcf out.indel.dp.q.recode.vcf --positions synnsyn.list --recode --out out.nsyn











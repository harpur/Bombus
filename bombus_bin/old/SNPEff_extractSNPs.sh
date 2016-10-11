
#Run SNPEFF for SNP information

java -jar /usr/local/lib/snpEff2/snpEff.jar Bimp -o txt /media/data1/bombus/vcf_Bimp/Bimp.indel.dp.bl.recode.vcf > Bimp.filt.snpeff.eff


sed '/WARNING/ !d' Bimp.filt.snpeff.eff > warnexons.filt.eff
sed '/WARNING/ d' Bimp.filt.snpeff.eff > Bimp.filt.1.snpeff.eff
sed '/INTERGENIC/ !d' Bimp.filt.1.snpeff.eff > interg.eff
sed '/UPSTREAM:/ !d' Bimp.filt.1.snpeff.eff > upstr.eff
sed '/INTRON/ !d' Bimp.filt.1.snpeff.eff > intron.eff
sed '/SYNONYMOUS/ !d' Bimp.filt.1.snpeff.eff > exons.filt.eff



#Run VCFTools for Depth Information

vcftools --vcf /media/data1/bombus/vcf_Bimp/Bimp.indel.dp.bl.recode.vcf --site-depth
# Bombus analysis archive



####contains all analyses performed on datasets in /bombus/


####Sample Information 
#####Bombus impatiens
Specimen ID|Collection ID|Collection location|Collection coordinates
------------ | ------------- | ------------- | -------------
Bimp01|08.06.13.01|G. Ross Lord Park|43°46.564' N, 79°27.994' W
Bimp02|08.06.13.02|Sunnybrook Park|43°43.379' N, 79°21.788' W
Bimp03|08.15.13.01|High park|43°38.702' N, 79°27.925' W
Bimp04|08.12.13.01|Boyd Apiary|43°49.503' N, 79°36.481' W
Bimp05|08.13.13.01|Rouge park|43°49.153' N, 79°10.234' W
Bimp06|08.14.13.01|Brickworks|43°41.183' N, 79°21.738' W
Bimp07|08.14.13.02|Eglington Park|43°42.332' N, 79°24.284' W
Bimp08|08.15.13.03|Eglington Flats|43°41.062' N, 79°30.221' W
Bimp09|08.16.13.01|Laurence's house|43°39.756' N, 79°26.263' W 
Bimp10|08.16.13.03|Lindylou park|43°44.897' N, 79°32.591' W 
Bimp02a|08.06.13.03|Sunnybrook Park|43°43.379' N, 79°21.788' W
Bimp09a|08.16.13.02|LP house|43°39.756' N, 79°26.263' W 
Bimp02b|08.06.13.03|Sunnybrook Park|43°43.379' N, 79°21.788' W
Bimp09b|08.16.13.01|Laurence's house|43°39.756' N, 79°26.263' W 
Bimp11|08.09.13.02|Rasberry farm|43°57.721' N, 79°33.270' W

####Alignment details
Alignment followed standard [GATK Best Practices for aligning genomes](http://www.broadinstitute.org/partnerships/education/broade/gatk-best-practices-and-building-analysis-pipelines-queue). In brief:
1) Trim fastq of adaptor sequences 
2) Align with BWA v0.7.5a-r405 sampe to the Bobmus terrestris reference genome (see below)
3) Remove Duplicate Reads with PCIARD v 1.141  AddOrReplaceReadGroups 

####SNP Calling
SNP calling was performed with GATK v 3.5-0-g36282e4 UnifiedGenotyper, as below, for all samples (see bombusalign.sh). We repeated this call with --glm indel for indels and again using --ploidy 2 to identify regions of the genome to exclude (see below)

<pre><code>gatk -R $BOMBUS/BimpFullGenome.fa -T UnifiedGenotyper \
	-I bams.list  \
	-o out.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-glm SNP  \
	-ploidy 1 &
</code></pre>

Following these raw SNP and indel calls, I trimmed all SNPs within 10bp of an indel and all SNPs with QD < 5.0, FS > 40.0, or MQ < 25.0. I also remove all SNPs within 5 bp of a SNP that was called as heterozygotic in at least 2 drone samples (putativeCNV.list)

####SNP Functional Roles
Identified putative functional role of SNPs using SNPEFF v XXXXXXXX. I collated a list of genes with polymorphic stop codons that will be masked from analysis





####GFF and Genome Versions 
We downloaded the GFF v 2.0 from NCBI and the corresponding reference genome. We use make_snpeff_input.py to take the longest transcript and if equal lengths, takes the first case.
Following SNPEff (below), we also excluded genes within comparisons that stop-codons or were of incorrect transcript length
I pulled all the CDS from the genome and translated them (with THIS SCRIPT). BLAST'ed them against mtDNA to ensure we didn't have any mtDNA sequence within the genome.  


####Masking the Bombus Genome
putativeCNV.list from heteroSNPs

mask via BLAST?







#####Conversions
For conversion between versions of the genome, orthologs, and gene names consult /GFF/ . We used BLAST best matches for conversion between genomes and for Fly "orthologs". For honey bee orthologs, we used a reciprocal best blast.
headers of the conversion file (BIMP_XP_GBoldnew_FBGN_recip.txt):
XP BIMP GBnew GBold FBGN description



####Admixture
I ran admixture on each species individually and then pooled together (see admixture.sh)
For Bimp K = 1 (CV error = 2.07322)
For Bterr K = 1
For Bmel K = 1


####Diversity
I ran BombusPi.r for each species, outputting diversity for each gene. 



####GO
Pep_DMELvsBIMPBlast.out
and Pep_DMELvsBIMPBlast.out




####UTR Gamma
Need to create a UTR GFF file. I did this from the original GFF file 















<!---
http://jeb.biologists.org/content/216/18/3474.full
http://rspb.royalsocietypublishing.org/content/281/1780/20132419.short
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131825/
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4074288/
http://www.pnas.org/content/108/18/7472.full




Quick info on the species used (From Williams...Colla et al. 2014):

B. melanopygus






regarding metabolis signatures:
https://books.google.ca/books?id=ifOcBAAAQBAJ&printsec=frontcover
http://onlinelibrary.wiley.com/doi/10.1111/mec.13410/epdf
http://www.ncbi.nlm.nih.gov/pubmed/26453894







-->
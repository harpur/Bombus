# Bombus analysis archive



####contains all analyses performed on datasets for Bombus Pop Geno (Harpur et al. 2016)

####Sample Information 
Sample|BioProject|Accession|Species|Location|Identified By
------------ | ------------- | ------------- | ------------- | ------------- | -------------
Bimp-10 |PRJNA347806|SAMN05894455|Bombus impatiens|Canada:Toronto:Lindylou park|Packer
Bimp-1  |PRJNA347806|SAMN05894456|Bombus impatiens|Canada:Toronto:G. Ross Lord Park|Packer
Bimp-11 |PRJNA347806|SAMN05894457|Bombus impatiens|Canada:Toronto:Rasberry farm|Packer
Bimp-2b |PRJNA347806|SAMN05894458|Bombus impatiens|Canada:Toronto:Sunnybrook Park|Packer
Bimp-3  |PRJNA347806|SAMN05894459|Bombus impatiens|Canada:Toronto:High park|Packer
Bimp-4|PRJNA347806|SAMN05894460|Bombus impatiens|Canada:Toronto:Boyd Apiary|Packer
Bimp-5|PRJNA347806|SAMN05894461|Bombus impatiens|Canada:Toronto:Rouge park|Packer
Bimp-6|PRJNA347806|SAMN05894462|Bombus impatiens|Canada:Toronto:Brickworks|Packer
Bimp-7|PRJNA347806|SAMN05894463|Bombus impatiens|Canada:Toronto:Eglington Park|Packer
Bimp-8|PRJNA347806|SAMN05894464|Bombus impatiens|Canada:Toronto:Eglington Flats|Packer
Bter-81  |PRJNA347806|SAMN05894465|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-82  |PRJNA347806|SAMN05894466|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-83  |PRJNA347806|SAMN05894467|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-84  |PRJNA347806|SAMN05894468|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-85  |PRJNA347806|SAMN05894469|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-87|PRJNA347806|SAMN05894470|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-88|PRJNA347806|SAMN05894471|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-89|PRJNA347806|SAMN05894472|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-90|PRJNA347806|SAMN05894473|Bombus terrestris|United Kingdom: Norwich|Bourke
Bter-91|PRJNA347806|SAMN05894474|Bombus terrestris|United Kingdom: Norwich|Bourke
BmelR003|PRJNA347806|SAMN05894475|Bombus melanopygus|USA:Oregon:Cape Sebastian State Park|Hines
R007|PRJNA347806|SAMN05894476|Bombus melanopygus|USA:Oregon:Rogue river, Jerry Flat Rd.|Hines
R001|PRJNA347806|SAMN05894477|Bombus melanopygus|USA:Oregon:Tahkenitch Creek Trailhead Rd|Hines


####Alignment details
Alignment followed previous standard [GATK Best Practices for aligning genomes](http://www.broadinstitute.org/partnerships/education/broade/gatk-best-practices-and-building-analysis-pipelines-queue). In brief:
Example details in bombusalign.sh

	1. Trim fastq of adaptor sequences 
	2. Align with BWA v0.7.5a-r405 sampe to the Bobmus terrestris reference genome (see below)
	3. Remove Duplicate Reads with PCIARD v 1.141  AddOrReplaceReadGroups 

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

Following these raw SNP and indel calls, I trimmed all SNPs within 10bp of an indel and all SNPs with QD < 5.0, FS > 40.0, or MQ < 25.0. I also remove all SNPs within 5 bp of a SNP that was called as heterozygotic in at least 2 drone samples (putativeCNV.list), and diallelic

####SNP Functional Roles
Identified putative functional role of SNPs using SNPEFF v 3.6c. [I created a custom B. impatiens database](http://snpeff.sourceforge.net/SnpEff_manual.html#databases) with the GFF and reference fasta as a text output (-o txt).

Final output for this is "exons.eff" for each species and "warnexons.eff" that lists genes with polymorphic stop codons. "exons.eff" was also used to compile a list of SNPs within genes for each species, called "out.nsyn.recode.vcf"


####GFF and Genome Versions 
We downloaded the GFF v 2.0 (/data/genes.gff for CDS regions) from NCBI and the corresponding reference genome. We use make_snpeff_input.py to take the longest transcript and if equal lengths, takes the first case.
Following SNPEff (below), we also excluded genes within comparisons that stop-codons or were of incorrect transcript length
I pulled all the CDS from the genome and translated them. BLAST'ed them against mtDNA to ensure we didn't have any mtDNA sequence within the genome.  


####Masking the Bombus Genome
1. BombusCNV.list contains a list of sites in the genome in which at least 2/10 haploid drone samples were called as heterozygotic in either of the Bombus genomes. Created from merged "putativeCNV.list" files.
2. PolyStop.list contains a list of genes with stop codons (from compiled "warnexons.eff") 
3. ScaffoldSummary contains statistics on scaffolds in BIMP genome 


####Running SnIPRE
[Bustamante's SnIPRE](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002806) was used to to estimate the selection coefficient on each gene inthe genome, along with MK test statistics. 

SnIPRE requires the following R packages and sources:

	1. lme4
	2. R2jags
	3. arm
	4. B_SnIPRE_source.R (from Bustamante)
	5. SnIPRE_source.R (from Bustamante)
	6. my.jags2.R (from Bustamante)
	7. SnIPRE.bug (from Bustamante)

	
The SnIPRE input file looks like the following table. This is an example, not real data, but headers must be the same

<pre><code>      GeneID GeneID PR FR PS FS Tsil Trepl nout npop
 NP_001267051.1  3 35  3 41 364.66667  1000   10    8
 NP_001267052.1  9 27  6 12 244.00000  1000   10    8
 XP_003484388.1 10 44  4  1  65.33333  1000   10    8
...
</code></pre>

Go get this output, I merged my vcf files for each species using vcf-merge:
<pre><code> 
vcf-merge out.nsyn.recode.vcf.gz out.nsyn.recode.vcf.gz | bgzip -c > /vcf/bombus.vcf.gz
</code></pre>

CreateMK.r builds an MK table from this merged VCF and it's corresponding SNPEFF txt output. It then compiles the MK table  into a table for SnIPRE(MK.snipre). As well, writes out MK test results and Neutrality index (MKresults). This test does not add +1 to every cell in the MK matrix to account for 0 entries. To Do: Make this script take arguments and run from command line.

The Neutrality index can be plotted as per [Li et al. (2008)](http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2008.00486.x/full) using MKplot.r (NI.pdf) 

SNIPRE is then run with SNIPRE.R
<pre><code> 
Rscript SNIPRE.R "MK.snipre"
</code></pre>

Bayesian output in .csv ending in ".bayesianresults" (in /data/)


####Data Analyses

#####Conversions
For conversion between versions of the genome, orthologs, and gene names consult /GFF/ . We used BLAST best matches for conversion between genomes and for Fly "orthologs" (blastp; e-value 1e-6). For honey bee "orthologs", we used a reciprocal best blast.
headers of the conversion file (BIMP_XP_GBoldnew_FBGN_recip.txt):
<pre><code> 
XP BIMP GBnew GBold FBGN description
</code></pre>


We also took advantage of [OrthoDB](http://orthodb.org/) and created the file "OrthoDB7_w_gamma.txt". This file outlines the likely taxonomic status of genes in the Bombus genome. Headers for the file are Genus, Bombus gene ID, taxonomic class, orthodb id, the number of copies (genera *2), the likely Apis ortholog, and that ortholog's gamma. 
<pre><code> 
genus	id	trg	og	n	ortho	gamma
</code></pre>

This file has only 1:1 orthologs, I've removed all more complex relationships. 

<!---
We also created the file "OrthoDB7_onetoone_Apis_Bombus.txt" that has the one-to-one orthologs from OrthoDB7 between Apis and Bombus and their respective gamma estimates. 

CHECK THAT!
-->


#####Woodard et al. 2014 data for Bombus DEGs
Hollis Woodard kindly provided the supplemental data from [her 2014 study](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4027384/) in which she compared expression (via microarray) between and among different Bombus castes in a very nice full factorial design. This study was made from Bombus ESTs (yay) and annotated against non-Bombus.  

I re-annotated by taking the probes and performing best BLASTN against all Bombus impatiens genes  (1e-6). I then inputted the BBM into Woodard et al's table (now called "rawDEGs.txt"). Some probes matched multiple genes. To deal with this, I took the mean-fold change in expression, mean gamma, and mean p-value from Woodard et al.'s data set for each gene with a duplicate. I repeated thsi aanalysis with only 1:1 matches (i.e. no duplicates) and got very similar results (see figure QWD_rep_Selection_NODUPLICATES.pdf). A gene was defined as differentially expressed in one caste/life history stage over another if it was significantly differentially expressed in all comparisons containing it but not sigificantly differentially expressed in any other comparisons. We defined a gene as ifferentially expressed in a reproductive caste if it was  in Foundresses or Queens and a gene is said to be differentially expressed in non-reproductive if it was differentially expressed in either workers or gynes. That is, we used the same definitions as [Woodard et al. (2014)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4027384/). The ANOVA model results from Woodard et al. can be found in the file "hollismodel". It looks like the example below where a 0 indicates not differentialyl expressed within that comparison, a 1 indicates over-expression in those comparisons and a -1 indicates udnerexpression in that comparison. For example, NP_001267052.1 is over expressed during brood care (MainBrood).

<pre><code> 
BimpGene	FW	QG	FQ	WG	MainRep	MainBrood	Interact
NP_001267052.1	0	0	0	1	0	1	0
NP_001267051.1	0	1	0	1	0	0	-1
XP_003484388.1	0	1	0	1	1	1	-1
XP_003484388.1	0	-1	0	-1	-1	-1	0
</code></pre>

 
 
#####Apis DEGs
For Apis mellifera DEGS I pulled the supplemental material from [Harpur et al. 2014](http://www.pnas.org/content/111/7/2614.abstract): "AMELDEGs.txt"




<!---

Analyses are stored in /data/ and not git.




####Admixture
I ran admixture on each species individually and then pooled together (see admixture.sh)
For Bimp K = 1 (CV error = 2.07322)
For Bterr K = 1
For Bmel K = 1

####GO
Pep_DMELvsBIMPBlast.out
and Pep_DMELvsBIMPBlast.out






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







####Fucking around with caste-spp genes:
I took my ortho list (1:1 amel to bombus)
I integrated the estimates of gamma in each spp
I also integrated the caste-specific gene expression patterns from Harpur et al. and Woodard et al. Comparing gene expression with protein expression seems silly...














-->
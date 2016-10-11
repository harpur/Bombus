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





####GFF
We downloaded the GFF v 2.0 from NCBI. We use make_snpeff_input.py to take the longest transcript and if equal lengths, takes the first case.
Following SNPEff (below), we also excluded genes within comparisons that stop-codons or were of incorrect transcipt length
I pulled all the CDS from the genome and translated them (with THIS SCRIPT). BLAST'ed them against mtDNA to ensure we didn't have any mtDNA sequence within the genome. None. 


#####Conversions
For conversion between versions of the genome, orthologs, and gene names consult /GFF/ . We used BLAST best matches for conversion between genomes and for Fly "orthologs". For honey bee orthologs, we used a reciprocal best blast.
headers of the conversion file (BIMP_XP_GBoldnew_FBGN_recip.txt):
XP BIMP GBnew GBold FBGN description

####Alignment details
1) This

2) This 




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




makeblastdb -in BIMP_CDS.fas -dbtype nucl -out BimpCDS -parse_seqids

#Nucl
 blastn -query bimp_OGSv1.0_cds.fa -db BimpCDS -out CDSvsBEEBASE.out  -outfmt "6 qseqid sseqid sstart send qstart qend length nident mismatch evalue" -evalue 1e-60

 
source("/media/data1/forty3/brock/scripts/ReferencetoDataframe.r")
dmel = Reftodf(data.frame="dmel-all-CDS-r6.04.fasta",nms="N")
chrom = gsub(".*;parent=","",dmel$chrom)
chrom = gsub(",FBtr.*","",chrom)
dmel$chrom = chrom
dmel=dmel[!duplicated(dmel$chrom),]

 

dmel = Reftodf(data.frame="Bimp_PEP.fas",nms="N")
dmel$seqs=gsub("-","",dmel$seqs)


makeblastdb -in out.fa -dbtype nucl -out DMELCDS -parse_seqids

 blastn -query Bimp_CDS.fas -db DMELCDS -out DMELvsBIMP.out  -outfmt "6 qseqid sseqid sstart send qstart qend length nident mismatch evalue" -evalue 1e-6
 
 
 #Protein
 
 cd /media/data1/forty3/brock/blast
 
 blastp -query Bimp_PEP.fas -db DMELpep -out DMELvsBIMP.out  -outfmt "6 qseqid sseqid sstart send qstart qend length nident mismatch evalue" -evalue 1e-6
 
 
best_blast.pl DMELvsBIMP.out
 
 
 
 
 
 
 #mtDNA blast
 cp /media/data1/bombus/BimpFullGenome.fa /media/data1/bombus/blast
 makeblastdb -in BimpFullGenome.fa -dbtype nucl -out BIMBPFULLGENOME -parse_seqids
 blastn -query mtDNABHYP.fas -db  BIMBPFULLGENOME -out DMELvsBIMP.out  -outfmt "6 qseqid sseqid sstart send qstart qend length nident mismatch evalue" 
 
 
 
 
 #mtDNA protein Blast
cd /media/data1/forty3/brock/blast
#sed -i 's/-//g' Bimp_PEP.fas 
makeblastdb -in Bimp_PEP.fas -dbtype prot -out BIMPPEP -parse_seqids
blastp -query BombusmtDNAprots.fas -db   BIMPPEP -out mtDNAvsBIMP.out  -outfmt "6 qseqid sseqid sstart send qstart qend length nident mismatch evalue" 

 
 
 
 
 
 
 
 
 
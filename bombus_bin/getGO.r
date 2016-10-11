rm(list=ls())
#library(ggplot2)
#library(gridExtra)
#library(ggbio)
#library(gtools)
#library(RMySQL)
library(GOstats)
library(GSEABase)

## ---------------- Function to calcualte GO terms and perform hypergeometric test againts background null distribution

get_GO <- function (infile){

   intr = read.delim(infile, sep="\t", header = F)
   genes = intr[intr$V6 == "gene",]
   rawgeneids = unique(droplevels(genes[genes$V6 == "gene",]$V12))
   geneids = sub("ID=", "", rawgeneids)

   outcomeBP <- hyperGTest(GSEAGOHyperGParams(name = "lrt", geneSetCollection=gsc,geneIds = unique(as.character(unlist(geneids, use.names=F))), universeGeneIds=as.character(universe),ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
   outcomeMF <- hyperGTest(GSEAGOHyperGParams(name = "lrt", geneSetCollection=gsc,geneIds = unique(as.character(unlist(geneids, use.names=F))), universeGeneIds=as.character(universe),ontology = "MF",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
   outcomeCC <- hyperGTest(GSEAGOHyperGParams(name = "lrt", geneSetCollection=gsc,geneIds = unique(as.character(unlist(geneids, use.names=F))), universeGeneIds=as.character(universe),ontology = "CC",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))

   towriteBP = cbind(infile, summary(outcomeBP))
   towriteMF = cbind(infile, summary(outcomeMF))
   towriteCC = cbind(infile, summary(outcomeCC))

   return(list("BP" = towriteBP, "MF" = towriteMF, "CC" = towriteCC))


}

## ----------------------- calcualte go temrs for significant snps -----------------------------

## read and process all go data
go = read.csv("go.csv")					# input file
universe <- unique(go$gene)
goFrame=GOFrame(go[go$gene %in% universe,],organism="Apis mellifera")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
GO=get_GO(infile = "significant.intr")    # input file
bp=GO$BP
mf=GO$MF
cc=GO$CC

## -------------------- permutation test at chromosome level ------------------------------------

## load all snps and significant snps
genome=read.delim("snps.bed", header=F)      # input file
genome$chr=sub("\\..*","",sub("Group","",genome$V1))
snps=read.delim("significant.bed", header=F)    # input file
snps$chr=sub("\\..*","",sub("Group","",snps$V1))

## frames to stores random snps and their go terms, and the number of permutations
rand.BP=c()
rand.MF=c()
rand.CC=c()
no.permut=1000

## sample for a limit times, and at each iteration select as much snps at each chromosome as orignially significant ones
for(i in 1:no.permut){

   frame=c()
   for(c in unique(snps$chr)){
      num=nrow(snps[snps$chr==c,])
      genome.sub=genome[genome$chr == c, ]
      rand=genome.sub[sample(nrow(genome.sub),num),]
      frame=rbind(frame, rand)
   }

   write.table(frame[,-4],"rand.bed",row.names=F, col.names=F, quote=F, append=F, sep="\t")
   system("intersectBed -a rand.bed -b amel_OGSv3.2.gff3 -wa -wb > rand.intr")

   GO.rand=get_GO(infile = "rand.intr")

   rand.BP=rbind(rand.BP,GO.rand$BP)
   rand.MF=rbind(rand.MF,GO.rand$MF)
   rand.CC=rbind(rand.CC,GO.rand$CC)

   print(i)

}

## calculate the frequency of each go term
freq.BP=table(rand.BP$GOBPID) / no.permut
freq.MF=table(rand.MF$GOMFID) / no.permut
freq.CC=table(rand.CC$GOCCID) / no.permut

## append the frequency of each origianl go term
bp=cbind(bp, "Pvalue.permutation"=freq.BP[bp$GOBPID])
mf=cbind(mf, "Pvalue.permutation"=freq.MF[mf$GOMFID])
cc=cbind(cc, "Pvalue.permutation"=freq.CC[cc$GOCCID])

## write the data to files
outb = "goBPterms.txt"
file.create(outb)

outm = "goMFterms.txt"
file.create(outm)

outc = "goCCterms.txt"
file.create(outc)

write.table(bp, file=outb, row.names=F, quote=F, append=F, sep=",", col.names=T)
write.table(mf, file=outm, row.names=F, quote=F, append=F, sep=",", col.names=T)
write.table(cc, file=outc, row.names=F, quote=F, append=F, sep=",", col.names=T)

## done
print("wrote go terms")

#!/usr/bin/Rscript

###
# Pi Across the Whole Genome in size W windows
###




#	Functions -------------------------------
source("/media/data1/forty3/brock/scripts/VCFFunctions.r")


# Load GFF,change to genes ---------------------------
gff = read.table(file = "/media/data1/bombus/GFF/BIMPPGFF_R.txt")
sts = aggregate(gff$V4, by=list(gff$V9), min)
ends = aggregate(gff$V5, by=list(gff$V9), max)
genes = gff[!duplicated(gff$V9),]
genes = genes[order(genes$V9),]
genes$V4 = sts$x
genes$V5 = ends$x

# Load VCF ---------------------------
args = commandArgs(trailingOnly = TRUE)
i = args[1]	
x =  Read.VCF()
#VCF files look like this 
#1.1	802	.	C	T	398.27	PASS	.	GT:AD:DP:GQ:PL	0/0:35,0:34:96:0,96,1243	0/0:19,0:19:57:0,57,743	0/0:26,0:25:75:0,75,967	0/0:26,0:25:75:0,75,951	0/0:30,0:29:75:0,75,985	0/0:39,0:38:99:0,111,1424	0/0:24,0:23:66:0,66,853	0/0:38,0:36:99:0,108,1412	0/0:30,0:29:84:0,84,1094	0/1:13,14:25:99:433,0,390	0/0:47,0:45:99:0,129,1652
#Alleles coded as 0/0, 0/1, 1/1


#	TD Pre-amble: -------------------------------
st=grep("FORMAT", unlist(x[1]))+1
x=x[-1]
en=length(unlist(x[1]))
#open df and get the number of samples (after FORMAT)
numInds=(en-st)+1
N=numInds
a1=0
a2=0
for(u in c(1:(N-1))){
a1 = c(a1,(1/u))
a2 = c(a2, 1/(u*u))
}
a1=sum(a1);a2=sum(a2)
b1 = (N+1)/(3*(N-1))#see Tajima 1989
b2 = (2*(N*N+N+3))/(9*N*(N-1))
c1 = b1 - 1/a1
c2 = b2-(N+2)/(a1*N)+a2/(a1*a1)
e1 = c1/a1
e2 = c2/(a1*a1+a2)

#	remove tri-allelic SNPs (x[5] is alternative allele) -------------------
len = sapply(x,function(x) nchar(x[5]))
x = x[len == 1]


#	Output SNP locations and scaffolds  -------------------
SNPs = sapply(x, function(x) rbind(x[1], x[2]))


for(i in 1:nrow(genes)){
#A for loop is really inefficient to do this, I'll vectorize later if I have to do it agian.
	#this gets out SNPs/exon
	snps = SNPs[2,][SNPs[1,] == genes[i,1]]
	
	#Gets all SNPs on the scaffold of the gene
	blah = outer(as.numeric(snps), as.numeric(genes[i,4]), ">=") 
	blah1 = outer(as.numeric(snps), as.numeric(genes[i,5]), "<=") 
	blah = (which(blah1 == "TRUE" & blah == "TRUE", arr.ind = T)) #SNPs is row, Gene is Col
	#37,47
	if(is.null(nrow(blah))){print("nope")}else{		
	if(nrow(blah)=="0"){print("nope")}else{	
		print(i)
		w = (genes[i,5]-genes[i,4]) #window size
		ends = 10 + (numInds-1) # element in snps that has that last individaul
		
# Count the number of genotypes for allelle freq ---------------------------		
		snps = x[SNPs[1,] == genes[i,1]]
		ones=sapply(snps[blah[,1]],function(x) (length(x[10:ends][grep("^1/1:", x[10:ends])]))) #homo ref
		twos=sapply(snps[blah[,1]],function(x) (length(x[10:ends][grep("^0/0:", x[10:ends])]))) #homo alt
		hets=sapply(snps[blah[,1]],function(x) (length(x[10:ends][grep("^0/1:", x[10:ends])]))) #homo alt
		miss=sapply(snps[blah[,1]],function(x) (length(x[10:ends][grep("[.]/[.]:", x[10:ends])]))) 
		miss=(abs(twos+ones+hets));miss=abs(miss-((en-st)+1))
		miss=((en-st)+1)-miss #now, gets sample size
		numInds=(en-st)+1
		
		blah=as.character(ones>=twos)
# Allelle freq ---------------------------			
		#If TRUE, 2*ones, if FALSE, 2*twos
		blah[which(blah=="TRUE")]	= (2 * twos[which(blah=="TRUE")] + hets[which(blah=="TRUE")])
		blah[which(blah=="FALSE")] =	(2 * ones[which(blah=="FALSE")] + hets[which(blah=="FALSE")])
	
		#blah that have miss==0 need to be removed.
		blah[miss=="0" | miss=="1"]=0 	
		blah=as.numeric(blah)
		blah=blah*(miss-(blah/2)) #N-J
		blah=round(blah/(miss*(miss-1)),4) #N(N-1)
# Calculate Pi and S ---------------------------			
		
		#need to divide this by length of exon.
		pi=round(sum(blah)/(w),8)

		#Geting Theta--Follow's AD's logic 
		S=rep("NA",length(ones))
		nS=rep("NA",length(ones))
		nS[miss==ones]="1";S[miss==ones]="0"
		nS[miss!=ones]="0";S[miss!=ones]="1"
		
	#Above, miss becomes the samplesize (after correcting for missing values)
	#If the site is fixed for a polymorphism, miss==ones or miss==twos, so I
	#add on to the nS vector a 1, allowing me to sum the number of sites as sum(nS) and
	#sum(S). Theta can then be calculated as S/(windowSize * a1)
	theta = Window.Theta(S=as.numeric(S),a1=a1,w=w)

# Tajima's D (1989) ---------------------------	
	TD=w*(pi) - w*(theta)
	TD[pi==0 | theta==0]=NA
	p1 = sum(as.numeric(S))*e1
	p2 = sum(as.numeric(S))*e2*(sum(as.numeric(S))-1)
	p1p2=sqrt(p1+p2)
	TD=TD/p1p2

	blah=cbind(as.character(genes[i,9]), pi, theta, TD, w )
	print(as.character(genes[i,9]))
	
	
	}
		write.table(blah,file=paste(args[1], ".GTD",sep=""),quote=F,row.names=F,col.names=F,append=T)
}
}	







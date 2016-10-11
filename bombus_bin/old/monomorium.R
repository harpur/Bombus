library(ggplot2)
library(RMySQL)
library(edgeR)
library(GOstats)
library(GSEABase)


mydb = dbConnect(MySQL(), user='monomorium', password='pharaonis', dbname='monomorium', host='ECOEVO.UNIT.OIST.JP')
contrasts <- dbGetQuery(mydb,"SELECT * FROM contrasts") # read a priori contrasts
contrasts <- contrasts[,sort(names(contrasts))]  # this makes the contrasts matrix compatible with design matrix, which is sorted by level
factors <- dbGetQuery(mydb,"SELECT * FROM factors")  #read treatments


	#     ____  ____________                     __           _     
	#    / __ \/ ____/ ____/  ____ _____  ____ _/ /_  _______(_)____
	#   / / / / / __/ __/    / __ `/ __ \/ __ `/ / / / / ___/ / ___/
	#  / /_/ / /_/ / /___   / /_/ / / / / /_/ / / /_/ (__  ) (__  ) 
	# /_____/\____/_____/   \__,_/_/ /_/\__,_/_/\__, /____/_/____/  
	#                                          /____/               

counts <- dbGetQuery(mydb,"SELECT isoform,MP1,MP2,MP3,MP4,MP5,MP6,MP7,MP8,MP9,MP10,MP11,MP12,MP13,MP14,MP15,MP16,MP17,MP18,MP19,MP20,MP21,MP22,MP23,MP24 FROM expected_counts")
row.names(counts) <- counts$isoform
counts <- subset(counts,select=-c(isoform))
fpkm <- dbGetQuery(mydb,"SELECT isoform,MP1,MP2,MP3,MP4,MP5,MP6,MP7,MP8,MP9,MP10,MP11,MP12,MP13,MP14,MP15,MP16,MP17,MP18,MP19,MP20,MP21,MP22,MP23,MP24  FROM fpkm")
row.names(fpkm) <- fpkm$isoform
fpkm <- subset(fpkm,select=-c(isoform))
keep=rowSums(fpkm>= 1)>= ncol(counts)/2  #select isoforms where at least half of the libraries have FPKM >1
#create design matrix, and label rows and columns according to our contrasts
design <- model.matrix(~0+factors$factor,data=counts)
rownames(design) <- colnames(counts)
colnames(design) <- names(contrasts)

dge <- DGEList(counts=round(counts[keep,]),group=factors$factor)   #apply filtering!
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge,design,verbose=TRUE)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge,design)

# plot NMDS as a sanity check

# mdsplot <- plotMDS(dge,top=500)
# plot(mdsplot, main="",xlab="Dimension 1",ylab="Dimension 2")
# text(mdsplot$cmdscale.out[,1],mdsplot$cmdscale.out[,2],factors$factor)

#testing contrasts
dge_list <- list() # this list stores results of likelihood ratio tests, and a list of differentially expressed genes, after multiple comparison correction
for (i in 1:nrow(contrasts)) {
	lrt <- glmLRT(fit,contrast=t(contrasts[i,]))
	dge_list[[i]] <- list(lrt = lrt, lrt_et = decideTestsDGE(lrt,p=0.05,adjust="BH"))
}

	#    __________     __                                         _      __                         __ 
	#   / ____/ __ \   / /____  _________ ___     ___  ____  _____(_)____/ /_  ____ ___  ___  ____  / /_
	#  / / __/ / / /  / __/ _ \/ ___/ __ `__ \   / _ \/ __ \/ ___/ / ___/ __ \/ __ `__ \/ _ \/ __ \/ __/
	# / /_/ / /_/ /  / /_/  __/ /  / / / / / /  /  __/ / / / /  / / /__/ / / / / / / / /  __/ / / / /_  
	# \____/\____/   \__/\___/_/  /_/ /_/ /_/   \___/_/ /_/_/  /_/\___/_/ /_/_/ /_/ /_/\___/_/ /_/\__/  
	                                                                                                  

#load go terms, and append 
go <- dbGetQuery(mydb,'SELECT * FROM blast2go WHERE go != ""')
go$evidence <- "ISS"
go <- go[,c("GO","evidence","isoform")]
universe <- unique(go$isoform)
goFrame=GOFrame(go[go$isoform %in% universe,],organism="Monomorium pharaonis")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())


#in the list of go terms "upregulated" mean there is more transcription in the contrast with the positive intex
go_list <- list()
for (i in 1:nrow(contrasts)) {
	go_list[[i]] <- list( 
		upregulated = hyperGTest(
			GSEAGOHyperGParams(name = paste("contrast ",i,"upregulated"),
				geneSetCollection=gsc,geneIds = intersect(rownames(dge_list[[i]]$lrt$table[dge_list[[i]]$lrt_et==1,]),universe),
				universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over")),
		downregulated = hyperGTest(
			GSEAGOHyperGParams(name = paste("contrast ",i,"downregulated"),
				geneSetCollection=gsc,geneIds = intersect(rownames(dge_list[[i]]$lrt$table[dge_list[[i]]$lrt_et==-1,]),universe),
				universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over")))
}

	#  _       __     _ __                                     __     __           _____ __   
	# | |     / /____(_) /____     ________  ____  ____  _____/ /_   / /_____     / __(_) /__ 
	# | | /| / / ___/ / __/ _ \   / ___/ _ \/ __ \/ __ \/ ___/ __/  / __/ __ \   / /_/ / / _ \
	# | |/ |/ / /  / / /_/  __/  / /  /  __/ /_/ / /_/ / /  / /_   / /_/ /_/ /  / __/ / /  __/
	# |__/|__/_/  /_/\__/\___/  /_/   \___/ .___/\____/_/   \__/   \__/\____/  /_/ /_/_/\___/ 
	#                                    /_/                                                  

#collect data from lists
dge_table <-  data.frame(et = numeric(0),contrast = numeric(0))   #indexes of dge genes
fc_table <-  data.frame()
goterms <-  matrix(nrow=0,ncol=ncol(summary(go_list[[1]]$upregulated))+2)
for (i in 1:nrow(contrasts)) {
	goterms <- rbind(goterms,data.frame(summary(go_list[[i]]$upregulated),contrast=i,type= "upregulated"))
	goterms <- rbind(goterms,data.frame(summary(go_list[[i]]$downregulated),contrast=i,type= "downregulated"))
	dge_table <- rbind(dge_table,data.frame(et=dge_list[[i]]$lrt_et,contrast=i))
	fc_table <- rbind(fc_table,data.frame(dge_list[[i]]$lrt$table,contrast=i))
}
write.csv(goterms,"go.csv")
write.csv(fc_table,"fc.csv")
write.csv(dge_table,"dge.csv")

###
# TRG Status of OrthoDB7  genes in Apis and Bombus
###


#Pulls in the orthoDB "dist" and "Group" databases and outputs, for Apis and Bombus, the TRG status of each gene AND their orthologs between both these groups:

# File: BOMBUSAPISOrthoDB7	
#		ID	class	og_id	BvAortholog
#		Bimp00001	INSECTA	INS015963	NA
#		Bimp00002	INSECTA	INS006491	GB44944
#		Bimp00003	INSECTA	INS023540	GB51735
#		Bimp00004	INSECTA	INS012111	GB50927
#		Bimp00005	INSECTA	INS005158	GB51119







# Conservation Status Analyses with ORTHODB CITES ---------------------------------
conv = read.table(file="/media/data1/bombus/OrthoDB6/Bee_Dist_Data.txt",header=F)
gens = read.table(file="/media/data1/bombus/OrthoDB6/Bee_Group_Data.txt",header=F)

names(conv)=c("node","og_id",
"species","genes",
"ACEPH",
"AECHI",
"AFLOR",
"AGAMB",
"AMELL",
"APISU",
"BIMPA",
"BMORI",
"BTERR",
"CCINC",
"CFLOR",
"DMELA",
"DNOVA",
"DPLEX",
"EMEXI",
"HLABO",
"HSALT",
"LALBI",
"LHUMI",
"MQUAD",
"MROTU",
"NVITR",
"PBARB",
"PHUMA",
"SINVI",
"TCAST")

#conv = conv[which(conv$BIMPA!="0"), ]
ortho = rep("NA", nrow(conv))
#conv = conv[-c(1:4)]


#BOMBUS and APIS Genes
ortho[rowSums(conv[which(names(conv) %in% c("BIMPA","BTERR", "AFLOR","AMELL"))]) > 1 ] = "BOMBUS"
#conv = conv[which(ortho=="NA"),]	
conv$BIMPA = NULL ; conv$BTERR = NULL
conv$AMELL = NULL ; conv$AFLOR = NULL

	
#APOIDEA Genes
ortho[rowSums(conv[which(names(conv) %in% c("AFLOR",
"AMELL",
"EMEXI",
"LALBI",
"MQUAD",
"MROTU",
"HLABO",
"DNOVA"))]) >= 1] = "APOIDEA"

									
#HYMENOPTERA						
conv$AFLOR=NULL;conv$AMELL=NULL;conv$EMEXI=NULL;conv$LALBI=NULL;conv$MQUAD=NULL;conv$MROTU=NULL
conv$HLABO=NULL;conv$DNOVA=NULL

ortho[rowSums(conv[which(names(conv) %in% c("ACEPH",
"AECHI",
"CCINC",
"CFLOR",
"HSALT",
"LHUMI",
"PBARB",
"SINVI",
"NVITR"))]) >= 1] = "HYMENOPTERA"


#HYMENOPTERA
conv$ACEPH=NULL
conv$AECHI=NULL
conv$CCINC=NULL
conv$CFLOR=NULL
conv$HSALT=NULL
conv$LHUMI=NULL
conv$PBARB=NULL
conv$SINVI=NULL
conv$NVITR=NULL
ortho[rowSums(conv[-c(1:4)]) >= 1] = "INSECTA"
conv$node = ortho







# Get out gene IDs -------------------------------
gens = gens[which(gens$V3=="BIMPA" | gens$V3=="AMELL" ),]
ag.gens = aggregate(gens$V4, by = list(gens$V4), length)
names(ag.gens) = c("V4","x")
gens = merge(gens, ag.gens, by = "V4")
names(gens)[3] = "og_id"
gens = merge(gens, conv, by = "og_id")
gens = gens[order(gens$V4),]

node.ord = gens$node
node.ord[node.ord=="APOIDEA"]=2
node.ord[node.ord=="BOMBUS"]=1
node.ord[node.ord=="HYMENOPTERA"]=3
node.ord[node.ord=="INSECTA"]=4
gens$node.ord = as.numeric(node.ord)

ag.gens = aggregate(gens$node.ord, by = list(gens$V4), max)
node.ord = ag.gens$x
node.ord[node.ord==2]="APOIDEA"
node.ord[node.ord==1]="BOMBUS"
node.ord[node.ord==3]="HYMENOPTERA"
node.ord[node.ord==4]="INSECTA"
node.ord[grep("GB", gen.list)] =  paste(node.ord[grep("GB", gen.list)], ".A",sep="")
node.ord[grep("Bimp*", gen.list)] =  paste(node.ord[grep("Bimp*", gen.list)], ".B",sep="")
node.ord[node.ord=="BOMBUS.B"] ="BOMBUS"
node.ord[node.ord=="BOMBUS.A"] ="APIS"
node.ord=gsub("[.].*","",node.ord)
ag.gens$x = node.ord




#identify BOMBUS and APIS orthologs:
gens.ins = gens[grep("INS*",gens$V2),]
gens.ins = gens.ins [gens.ins$V3 %in% c("BIMPA", "AMELL"),]

write.list(gens.ins, file="/media/data1/bombus/workingFiles/ORTHODB6")
write.list(ag.gens, file="/media/data1/bombus/workingFiles/ORTHODB6_BIMP")


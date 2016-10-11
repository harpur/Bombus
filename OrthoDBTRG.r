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


#
cat ODB8* > ODBNEW
cut -f 1,2,3,4,5 ODBNEW > ODB
cut -f 1,5 ODB > spplist
#




# Load in Data frames -------------------------------
spp = read.table(file="spplist",header=T,sep = "\t")
odb = read.table(file="ODB", sep = "\t", header=T)

# Clean and quantify Species list -----------------------------
spp = spp[which(spp$organism != "organism"),] 
lvl = as.character(spp$odb8_level)
lvl[lvl == "34735:Apoidea"] = 1
lvl[lvl == "50557:Insecta"] = 3
lvl[lvl == "7399:Hymenoptera"] = 2
spp$lvl = as.numeric(lvl)
spp.list = aggregate(spp$lvl, by = list(spp$organism), min)
names(spp.list)[1] = "organism"

# Classify Genes -----------------------
odb = merge(odb, spp.list, by = "organism")
min.class = aggregate(odb$x, by = list(odb$odb8_og_id), max)
names(min.class) = c("odb8_og_id","trg")
odb = merge(odb, min.class, by = "odb8_og_id")

# Classify BIMP Genes -----------------------
bimp.odb = odb[grep("Bimp*",odb$protein_id),]
names(bimp.odb)[4] = "GB"

amel.odb = odb[grep("GB*",odb$protein_id),]
names(amel.odb)[4] = "GB"
amel.odb$GB = gsub("-PA","", amel.odb$GB)






#pull in bomb.gamma and amel.gamma
amel.gamma = read.table(file="/data2/bombus/workingFiles/gamma",header=T) #from Harpur PNAS 2014
amel.gamma = merge(amel.gamma, amel.odb, by = "GB")

amel.gamma1 =  aggregate(amel.gamma$trg, by = list(amel.gamma$GB), max)
amel.gamma1$gamma = aggregate(amel.gamma$gamma, by = list(amel.gamma$GB), mean)$x






conv = read.table(file = "/data2/bombus/git/data/BIMP_XP_GBoldnew_FBGN_recip.txt")
names(conv)[2] = "GB"

bimp.odb = merge(conv, bimp.odb, by = "GB")
names(bimp.odb)[c(1,2)] = c("BID","GB")



bomb.gamma = merge(gm.bter, gm.bmel, by="GeneID",suffix=c(".bter",".bmel"), all=T)
bomb.gamma = bomb.gamma[c(1,21,45)]
bomb.gamma = bomb.gamma[complete.cases(bomb.gamma),]
bomb.gamma = bomb.gamma[c(1,2)]
names(bomb.gamma) = c("GB", "gamma")
bomb.gamma = merge(bomb.gamma, bimp.odb, by = "GB", all.x=T)
bomb.gamma = bomb.gamma[complete.cases(bomb.gamma),]

bomb.gamma1 =  aggregate(bomb.gamma$trg, by = list(bomb.gamma$GB), max)
bomb.gamma1$gamma = aggregate(bomb.gamma$gamma, by = list(bomb.gamma$GB), mean)$x

lvl = as.character(bomb.gamma1$x)
lvl[lvl == "1"] = "Apoidea"
lvl[lvl == "3"] = "Insecta"
lvl[lvl == "2"] = "Hymenoptera"
bomb.gamma1$lvl = lvl










#pull in bomb.gamma and amel.gamma
conv = read.table(file = "/data2/bombus/git/data/BIMP_XP_GBoldnew_FBGN_recip.txt")
names(conv)[2] = "GB"

bimp.odb = merge(conv, bimp.odb, by = "GB")
names(bimp.odb)[c(1,2)] = c("BID","GB")



bomb.gamma = merge(gm.bter, gm.bmel, by="GeneID",suffix=c(".bter",".bmel"), all=T)
bomb.gamma = bomb.gamma[c(1,21,45)]
bomb.gamma = bomb.gamma[complete.cases(bomb.gamma),]
bomb.gamma = bomb.gamma[c(1,2)]
names(bomb.gamma) = c("GB", "gamma")
bomb.gamma = merge(bomb.gamma, bimp.odb, by = "GB", all.x=T)
bomb.gamma = bomb.gamma[complete.cases(bomb.gamma),]

bomb.gamma1 =  aggregate(bomb.gamma$trg, by = list(bomb.gamma$GB), max)
bomb.gamma1$gamma = aggregate(bomb.gamma$gamma, by = list(bomb.gamma$GB), mean)$x

lvl = as.character(bomb.gamma1$x)
lvl[lvl == "1"] = "Apoidea"
lvl[lvl == "3"] = "Insecta"
lvl[lvl == "2"] = "Hymenoptera"
bomb.gamma1$lvl = lvl





















#create barplot df
gamma.df1= bomb.gamma1#[which(bomb.gamma1$gamma>0),]
melted = melt(gamma.df1, id.vars=c("lvl", "Group.1"),measure.vars=c("gamma"))
means.sem = ddply(melted, c("lvl","variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)



means.sem = means.sem[which(means.sem$x!="Genera"),]
means.sem = means.sem[complete.cases(means.sem)]
means.sem

table(bomb.gamma1$lvl[bomb.gamma1$gamma>1]);table(bomb.gamma1$lvl)






















# Conservation Status Analyses with ORTHODB CITES ---------------------------------
conv = read.table(file="/data2/bombus/OrthoDB/Bee_Dist_Data.txt",header=F)
gens = read.table(file="/data2/bombus/OrthoDB/Bee_Group_Data.txt",header=F)

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
names(gens) =c("node","og_id","sp","GID","V5")

conv.save = conv
conv = merge(conv, gens, by = "og_id")
conv = conv[which(conv$sp=="AMELL" ),] #| conv$sp=="AMELL" 

#nrow(conv[conv$AMELL<2 & conv$BIMPA<2,]) #most restrictive case, single copy orthos
#conv = conv[which(conv$AMELL<2 & conv$BIMPA<2),] #most restrictive case, single copy orthos


#conv = conv[which(conv$BIMPA!="0"), ]
ortho = rep("NA", nrow(conv))
#conv = conv[-c(1:4)]


#BOMBUS and APIS Genes
ortho[rowSums(conv[which(names(conv) %in% c("AFLOR","AMELL"))]) > 1 ] = "SPP" #, "AFLOR","AMELL", "BIMPA","BTERR"
#conv = conv[which(ortho=="NA"),]	
#conv$BIMPA = NULL ; conv$BTERR = NULL
#conv$AMELL = NULL ; conv$AFLOR = NULL

	
#APOIDEA Genes
ortho[rowSums(conv[which(names(conv) %in% c("BIMPA","BTERR", ##OR....above
"EMEXI",
"LALBI",
"MQUAD",
"MROTU",
"HLABO",
"DNOVA"))]) >= 1] = "APOIDEA"

									
#HYMENOPTERA						
#conv$AFLOR=NULL;conv$AMELL=NULL;conv$EMEXI=NULL;conv$LALBI=NULL;conv$MQUAD=NULL;conv$MROTU=NULL
#conv$HLABO=NULL;conv$DNOVA=NULL

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
#conv$ACEPH=NULL
#conv$AECHI=NULL
#conv$CCINC=NULL
#conv$CFLOR=NULL
#conv$HSALT=NULL
#conv$LHUMI=NULL
#conv$PBARB=NULL
#conv$SINVI=NULL
#conv$NVITR=NULL

ortho[rowSums(conv[which(names(conv) %in% c("AGAMB", "APISU", "BMORI", "DMELA", "DPLEX", "PHUMA", "TCAST"))]) >= 1] = "INSECTA"
conv$node = ortho


node.ord = conv$node
node.ord[node.ord=="APOIDEA"]=2
node.ord[node.ord=="SPP"]=1
node.ord[node.ord=="HYMENOPTERA"]=3
node.ord[node.ord=="INSECTA"]=4
conv$node.y = as.numeric(node.ord)



ag.gens = aggregate(conv$node.y, by = list(conv$GID), max)
ag.gens$len = aggregate(conv$node.y, by = list(conv$GID), length)$x
node.ord = ag.gens$x
node.ord[node.ord==2]="Apoidea"
node.ord[node.ord==1]="Genera"
node.ord[node.ord==3]="Hymenoptera"
node.ord[node.ord==4]="Insecta"
ag.gens$x = node.ord

apis.trg = ag.gens # and bimp.trg


names(apis.trg)[1] = "GB"; apis.trg = apis.trg[-3]
names(bimp.trg)[1] = "GB"






#pull in bomb.gamma and amel.gamma
conv = read.table(file = "/data2/bombus/git/data/BIMP_XP_GBoldnew_FBGN_recip.txt")
names(conv)[2] = "GB"
bimp.trg = merge(conv, bimp.trg, by = "GB", all.x = T)
bimp.trg = bimp.trg[c(2,7)]; names(bimp.trg)[1] = "GB"
ag.gens = rbind(bimp.trg, apis.trg)


bomb.gamma = merge(gm.bter, gm.bmel, by="GeneID",suffix=c(".bter",".bmel"), all=T)
bomb.gamma = bomb.gamma[c(1,21,45)]
bomb.gamma = bomb.gamma[complete.cases(bomb.gamma),]
bomb.gamma = bomb.gamma[c(1,2)]
names(bomb.gamma) = c("GB", "gamma")
bomb.gamma$genus = rep("Bumble bees", nrow(bomb.gamma))
amel.gamma$genus = rep("Honey bees", nrow(amel.gamma))
gamma.df = rbind(bomb.gamma, amel.gamma)
gamma.df = merge(gamma.df, ag.gens, by = "GB", all.x=T)
gamma.df = gamma.df[!duplicated(gamma.df$GB),]


#create barplot df

gamma.df1= gamma.df[which(gamma.df$gamma>0),]
melted = melt(gamma.df1, id.vars=c("genus","x", "GB"),measure.vars=c("gamma"))
means.sem = ddply(melted, c("genus","x","variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem = means.sem[which(means.sem$x!="Genera"),]
means.sem = means.sem[complete.cases(means.sem)]
means.sem























#plot

df = means.sem
#means.sem = df[df$n=="2",] #2 supports 

p<-ggplot(df, aes(x = factor(x), y = mean,  fill = factor(x))) + 
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	facet_grid(. ~ genus) +
	geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 	+
	scale_fill_manual(values=wes_palette(n=3, name="Chevalier")) +
	theme_classic() + 
	theme(
		axis.line.x = element_line(size=1.1),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_text(size = 14),
		axis.title.y = element_text(size = 14),
		axis.text.x=element_blank(),
		axis.text=element_text(size=13),
		axis.ticks = element_blank(),
		legend.position= c(0.2,0.9),
		strip.background = element_rect(colour="white", fill="white")
		) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	guides(fill=guide_legend(title=NULL)) +
	labs(
	x= "",
	y = expression(gamma)) 
p


































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
node.ord[node.ord=="SPP"]=1
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


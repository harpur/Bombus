###
# Bombus Comparative Genomics Figures 
###



# Load Functions and Dataframes --------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(wesanderson)
#http://www.ucl.ac.uk/~zctpep9/Archived%20webpages/Cookbook%20for%20R%20%C2%BB%20Colors%20(ggplot2).htm
#http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
source(file="/data2/bombus/git/GenomeR/VarFunct.r")
gm.bter = read.csv("/data2/bombus/data/MK.snipre.bayesianresults", header=T)
amel.gamma = read.table(file="gamma",header=T) #from Harpur PNAS 2014
#names(gm.bter)[21] = "gm.bter"
gm.bmel = read.csv("/data2/bombus/data/MKresults.bmel.bayesianresults", header=T)
conv = read.table(file="BIMP_XP_GBoldnew_FBGN_recip.txt",header=F);names(conv)[c(1,3)] = c("XP","GB")
ortho = read.table("OrthoDB7_onetoone_Apis_Bombus.txt",header=T)
ortho = ortho[which(complete.cases(ortho$b.gamma)),]
ortho = ortho[which(complete.cases(ortho$gamma)),]
degs = read.table(file="/data2/bombus/git/data/rawDEGs.txt",header=T) 
amel = read.table("/data2/bombus/git/data/AMELDEGs.txt", header=T)
ortho.gamma = read.table(file="OrthoDB7_w_gamma.txt",header=T,sep="\t")	
hol.model = read.table(file="hollismodel",header=T); names(hol.model)[1]="XP"
mk = read.table(file="MKresults",header=T); names(hol.model)[1]="XP"
conv = read.table(file="/data2/bombus/GFF/BIMP_XP_GBoldnew_FBGN_recip.txt",header=F)

degs = read.table(file="/data2/bombus/DEG/rawDEGs.txt",header=T) #this is now overwritten with gamma
amel = read.table("/data2/bombus/workingFiles/AMELDEGs.txt", header=T) #load in sd5 from PNAS


# Build gamma DF -----------------------------------------		
	#Here, I'm only taking genes with estimates in all Bombus species (N=9983)
bomb.gamma = merge(gm.bter, gm.bmel, by="GeneID",suffix=c(".bter",".bmel"), all=T)
bomb.gamma = bomb.gamma[c(1,21,45)]
bomb.gamma = bomb.gamma[complete.cases(bomb.gamma),]
bomb.gamma = bomb.gamma[c(1,2)]




#Histogram Gamma in Apis and Bombus --------------------------
# Build Histogram DF -----------------------------------------	
test = bomb.gamma; names(test) = c("name","gamma")
test$genus = rep("Bumble bees", nrow(test))
names(amel.gamma) = c("name","gamma")
amel.gamma$genus = rep("Honey bees", nrow(amel.gamma))
gamma.hist = rbind(amel.gamma, test)


#Plot comparative gamma Hist ---------------------------------	
pdf(file="Figure1_ApisBombusGammaHistogram.pdf",width=16,height=10)
p<-ggplot(gamma.hist, aes(x=gamma))+
	geom_histogram(color="black", fill="white", binwidth = 0.2)+
	facet_grid(genus ~ .) +
	theme_tufte() + 
	theme(
		strip.text.y = element_text(size = 16),
		axis.title.x = element_text(size = 16),
		axis.text=element_text(size=10)	) +
	labs(
		y= "",
		x = expression(gamma))
p
dev.off()

png("Figure1_ApisBombusGammaHistogram.png")
p
dev.off()




#Apis vs Bombus Ortholog Gamma (based on BBM) --------------------------
# Load data frame
#pdf(file="ApisvsBombusBBMGamma.pdf")

test = merge(gm.bter, gm.bmel, by="GeneID",suffix=c(".bter",".bmel")) 
test = test[c(1,21,10)]
names(test)[1] = names(conv)[1] = "XP"
test = merge(test, conv, by = "XP")
test = test[complete.cases(test),]
names(test)[5]="GB";names(amel.gamma)[1] = "GB"
test = merge(test, amel.gamma, by = "GB")

pdf(file="Figure2_ApisvsBombusOrthoDBGamma.pdf")
p<-ggplot(test, aes(x=BSnIPRE.gamma.bter, y=gamma)) +
	geom_point(shape=19,      # Use solid circles
			   alpha=1/2) + 	  # Make transparent
	coord_cartesian(ylim = c(-2, 11), xlim = c(-2, 6)) + 
	theme_few() +	
	theme(
			axis.text.y=element_text(size=12),
			axis.text.x=element_text(size=12)		
			
		) +
	labs(x = expression("Bumble bees" ), y = expression("Honey bees")) + 
	geom_abline(slope = 1, lty=2, col ="red")
p
dev.off()

png("Figure2_ApisvsBombusOrthoDBGamma.png")
p
dev.off()




#Apis vs Bombus Ortholog Gamma (Based on OrthoDB)--------------------------
# Load data frame
#ortho = read.table("/data2/bombus/workingFiles/OrthoGammaBIMPvAPIS.txt",header=T)
ortho = read.table("/data2/bombus/workingFiles/OrthoDB7_onetoone_Apis_Bombus.txt",header=T,sep="\t")


ortho = ortho[which(complete.cases(ortho$b.gamma)),]
ortho = ortho[which(complete.cases(ortho$gamma)),]



p<-ggplot(ortho, aes(x=b.gamma, y=gamma)) +
	geom_point(shape=19,      # Use solid circles
			   alpha=1/2) + 	  # Make transparent
	coord_cartesian(ylim = c(-2, 11), xlim = c(-2, 6)) + 
	theme_classic() +	
	theme(
			axis.text.y=element_text(size=10),
			axis.text.x=element_text(size=10)		
			
		) +
	labs(x = expression(italic("Bombus ")~gamma), y = expression(italic("Apis ")~gamma)) + 
	geom_abline(xintercept = -2, lty=2, col ="red")
p


pdf(file="ApisvsBombusOrthoDBGamma.pdf")
p
dev.off()


png(file="ApisvsBombusOrthoDBGamma.pdf")
p
dev.off()




#Gamma in Apis vs Bombus QWD genes --------------------------
#load dataframe
amel = merge(amel, amel.gamma, by ="GB")
#create labels	
degs = aggregate(degs[c(2:10)], by=list(degs$XP), mean)
names(degs)[1]="XP"
sp = rep("NDEG", nrow(degs))
sp[which(degs$FW<0.05 & degs$QG>0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "F" #F
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "Q" #Q
sp[which(degs$FW<0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG<0.05)] = "W" #W
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ>0.05 & degs$WG<0.05)] = "G"#G
degs$sp=sp
sp = as.character(degs$sp)
sp[sp=="F"] = "Reproductive"
sp[sp=="Q"] = "Reproductive"
sp[sp=="W"] = "Non-reproductive" #non-reproductive
sp[sp=="G"] = "Non-reproductive"
degs$sp2 = sp
sp = as.character(amel$Caste)
sp[sp=="Q"] = "Reproductive"
sp[sp=="W"] = "Non-reproductive"
amel$sp2 = sp
degs.plot = degs[c(1,2,12)]; degs.plot$genus = rep("Bumble bees", nrow(degs)); names(degs.plot)[1]="GB"
amel.plot = amel[c(1,3,4)]; amel.plot$genus = rep("Honey bees", nrow(amel.plot))

qwd.plot = rbind(amel.plot, degs.plot)
qwd.plot = qwd.plot[complete.cases(qwd.plot ),] #make this the df that I load.

#create barplot df
melted = melt(qwd.plot, id.vars=c("genus","sp2"),measure.vars=c("gamma"))
means.sem = ddply(melted, c("genus","sp2", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)

ann.text =  means.sem[c(1,2,4)]
ann.text$mean = c("18.6%","20.7%","31.0%","14.4%","29.8%","17.8%") #proportion >1 
ann.text$x = c(1,2,3, 1,2,3)
ann.text$y = c(0.3,0.3,0.3,0.15,0.3,0.15)


#plot 

p<-ggplot(means.sem, aes(x = factor(sp2), y = mean,  fill = factor(sp2))) +  
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, .9)) + 
	facet_grid(. ~ genus) +
	geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 	+
	scale_fill_manual(values=wes_palette(n=3, name="Chevalier")) +
	
	geom_text(data = ann.text,aes(x = x, y = y, label = mean)) +
	
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

pdf(file="Figure4_QWD_rep_Selection.pdf")
p
dev.off()
#http://stackoverflow.com/questions/24578352/add-a-segment-only-to-one-facet-using-ggplot2


ann_text <- data.frame(mpg = 15,wt = 5,lab = "Text",
                       cyl = factor(8,levels = c("4","6","8")))
p + geom_text(data = ann_text,label = "Text")





png(file="Figure4_QWD_rep_Selection.png")
p
dev.off()












# Orthology and TRG Status vs Gamma ------------------------------
#load dataframe
ortho = read.table(file="/data2/bombus/git/data/OrthoDB7_w_gamma.txt",header=T,sep="\t")	
	#list of all genes with direct orthologs between Apis and Bombus, based on orthoDB7. Eliminated all more complex relationships. 
gen = ortho$genus



#create barplot df
melted = melt(ortho, id.vars=c("genus","trg","n"),measure.vars=c("gamma"))
means.sem = ddply(melted, c("genus","trg", "n","variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)



#plot

df = means.sem
means.sem = df[df$n=="2",] #2 supports 

p<-ggplot(means.sem, aes(x = factor(trg), y = mean,  fill = factor(trg))) + 
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, 1.5)) + 
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





pdf(file="Figure3_TRGOrthologs.pdf")
p
dev.off()


png(file="Figure3_TRGOrthologs.png")
p
dev.off()



















#Supplemental Figs --------------------------------------------------------------------------------

			# Gamma between Orthologs (OrthoDB) ------------------
			# Load data 
			#ortho = read.table(file="/media/data1/bombus/workingFiles/OrthoGammaBIMPvAPIS.txt",header=T)
			#
			#	#plot	
			#	pdf(file="ApisvsBombusOrthoDBGamma.pdf")
			#	ggplot(ortho, aes(x=b.gamma, y=ortho.gamma)) +
			#		geom_point(shape=19,      # Use solid circles
			#				   alpha=1/5) + 	  # Make transparent
			#	labs(x = expression(italic("Bombus ")~gamma), y = expression(italic("Apis ")~gamma)) +
			#	theme_bw() +	
			#	theme(
			#		axis.text=element_text(size=12),
			#		#axis.ticks.x = element_blank(),
			#		axis.title=element_text(size=16,face="bold"),
			#		#axis.text.x = element_blank(),	
			#		panel.border = element_rect(size=1)
			#		) 
			#
			#	dev.off()




# Reproductive vs Brood Care Genes (Hollis' Paper) ----------------
#Load data
hol.model = read.table(file="/data2/bombus/DEG/hollismodel.txt",header=T)
test =merge(degs, hol.model, by="XP")
test = test[!duplicated(test$XP),]#choose randomly


#create barplot df
melted = melt(test, id.vars="MainRep",measure.vars=c("gamma"))
means.sem = ddply(melted, c("MainRep", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem$variable = rep("Reproductive")
names(means.sem)[1]="exp"

melted = melt(test, id.vars="MainBrood",measure.vars=c("gamma"))
means2.sem = ddply(melted, c("MainBrood", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means2.sem = transform(means2.sem, lower=mean-sem, upper=mean+sem)
means2.sem$variable = rep("Brood")
names(means2.sem)[1]="exp"
means.sem = rbind(means.sem,means2.sem)
#exp = means.sem$exp 
#exp[exp=="-1"] = "Under-expressed"
#exp[exp=="0"] = "NDEG"
#exp[exp=="1"] = "Over-Expressed"
#means.sem$exp = exp


#plot 
pdf(file="ReproductiveGenesBarPlot.pdf")



p<-ggplot(means.sem, aes(x = factor(exp), y = mean,  fill = factor(exp))) + 
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, .9)) + 
	facet_grid(. ~ variable) +
	geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 	+
	theme_classic() + 
	theme(
		strip.text.x = element_text(size = 12.5),
		axis.title.y = element_text(size = 12.5),
		axis.text.x=element_blank(),
		axis.text=element_text(size=10),
		axis.ticks = element_blank(),
		legend.position= c(0.2,0.9)
		) +
	
	scale_fill_manual(
		name = "",
		values = brewer.pal(n = 3, name = "Blues"),
		breaks=c("-1", "0", "1"),
        labels = c("Under-expressed",
                   "NDEG",
                   "Over-expressed")) +
	labs(
	x= "",
	y = expression(gamma)) 
	
p


dev.off()












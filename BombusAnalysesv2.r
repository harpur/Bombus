###
# Bombus Analyses
###








#NOTE: all in /media/data1/bombus until BOMBUS server updated. THis will be moved to /data







# Load Functions and Dataframes --------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
#http://www.ucl.ac.uk/~zctpep9/Archived%20webpages/Cookbook%20for%20R%20%C2%BB%20Colors%20(ggplot2).htm
#http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
source(file="/data2/bombus/git/GenomeR/VarFunct.r")
gm.bter = read.csv("MK.snipre.bayesianresults", header=T)
amel.gamma = read.table(file="gamma",header=T) #from Harpur PNAS 2014
	#names(gm.bter)[21] = "gm.bter"
	#gm.bmel = read.csv("/media/data1/bombus/Snipre/snipre_run/bayesian_results_BterToBimp.csv", header=T)
conv = read.table(file="BIMP_XP_GBoldnew_FBGN_recip.txt",header=F);names(conv)[c(1,3)] = c("XP","GB")
ortho = read.table("OrthoDB7_onetoone_Apis_Bombus.txt",header=T)
ortho = ortho[which(complete.cases(ortho$b.gamma)),]
ortho = ortho[which(complete.cases(ortho$gamma)),]
degs = read.table(file="rawDEGs.txt",header=T) 
amel = read.table("AMELDEGs.txt", header=T)
ortho.gamma = read.table(file="OrthoDB7_w_gamma.txt",header=T,sep="\t")	


	
					
				##################	
				#Exploration and confirmation
				table(gm.bter$BSnIPRE.class)
				table(gm.bter$BSnIPRE.Rclass)
				gm.bter$class = as.factor(paste(gm.bter$BSnIPRE.class, gm.bter$BSnIPRE.Rclass,sep=""))
				plot(gm.bter$BSnIPRE.Rest, gm.bter$BSnIPRE.gamma, col = gm.bter$class)

				#
				#names(gm.bter)[21] = "gm.bmel"
				test = merge(gm.bter, gm.bmel, by="Group",suffix=c(".bter",".bmel"))
				#test = test[c(1,21,45)]
				#Gamma correlated between species
				with(test,cor.test(BSnIPRE.gamma.bter,BSnIPRE.gamma.bmel))

				#        Pearson's product-moment correlation
				#
				#data:  BSnIPRE.gamma.bter and BSnIPRE.gamma.bmel
				#t = 85.815, df = 10006, p-value < 2.2e-16
				#alternative hypothesis: true correlation is not equal to 0
				#95 percent confidence interval:
				# 0.6396889 0.6622644
				#sample estimates:
				#      cor
				#0.6511206
				#
				#GC3?
				gc3 = read.table(file="/media/data1/bombus/workingFiles/ThirdcodonACGT.txt",header=T)
				test = merge(test, gc3, by="Group")
				plot(test$BSnIPRE.gamma.bter , test$GC3)	
				cor.test(test$BSnIPRE.gamma.bter , test$GC3)

				plot(test$BSnIPRE.gamma.bter[test$BSnIPRE.gamma.bter>1] , test$GC3[test$BSnIPRE.gamma.bter>1])	
				cor.test(test$BSnIPRE.gamma.bter[test$BSnIPRE.gamma.bter>1] , test$GC3[test$BSnIPRE.gamma.bter>1])	
				cor.test(test$BSnIPRE.gamma.bter[test$BSnIPRE.class.bter=="pos"] , test$GC3[test$BSnIPRE.class.bter=="pos"])


				#which genes have highest residuals?
				gm.lm = lm(BSnIPRE.gamma.bter~BSnIPRE.gamma.bmel, data=test) 
				gm.res = resid(gm.lm)
					#interesting, no gene with high residuals has oposite ev.pattern.

					



				test = gm.bter[c(1,21)]; names(test) = c("name","gamma")
				test$genus = rep("Bumblebees", nrow(test))
				names(amel.gamma) = c("name","gamma")
				amel.gamma$genus = rep("Honey bees", nrow(amel.gamma))	
					
				nrow(amel.gamma); nrow(amel.gamma[which(amel.gamma$gamma>1),])
				nrow(test); nrow(test[which(test$gamma>1),])	




				###################



#Histogram Gamma in Apis and Bombus --------------------------



# Build Histogram DF -----------------------------------------	
test = gm.bter[c(1,21)]; names(test) = c("name","gamma")
test$genus = rep("Bumblebees", nrow(test))
names(amel.gamma) = c("name","gamma")
amel.gamma$genus = rep("Honey bees", nrow(amel.gamma))
gamma.hist = rbind(amel.gamma, test)


#Plot comparative gamma Hist ---------------------------------	
pdf(file="ApisBombusGammaHistogram.pdf",width=16,height=10)
p<-ggplot(gamma.hist, aes(x=gamma))+
	geom_histogram(color="black", fill="white", binwidth = 0.2)+
	facet_grid(genus ~ .) +
	theme_few() + 
	theme(
		strip.text.y = element_text(size = 12),
		axis.title.x = element_text(size = 12),
		axis.text=element_text(size=10)	) +
	labs(
		y= "",
		x = expression(gamma))
p
dev.off()




#Gamma in Apis vs Bombus Orthologs (based on BBM) --------------------------

# Build plot DF -----------------------------------------	
test = gm.bter[c(1,21)]; names(test) = c("XP","BSnIPRE.gamma.bter")
names(amel.gamma) = c("GB","gamma")

	#test = merge(gm.bter, gm.bmel, by="Group",suffix=c(".bter",".bmel")) 
	#test = test[c(1,21,10)]
	#test = test[c(1,21)]
test = merge(test, conv, by = "XP")
test = test[complete.cases(test),]
names(test)[4]="GB";
#amel.gamma = read.table(file="/data2/bombus/git/bomb/gamma",header=T)
test = merge(test, amel.gamma, by = "GB")


#Plot Apis gamma vs Bombus ---------------------------------	
pdf(file="ApisvsBombusOrthoDBGamma.pdf")
p<-ggplot(test, aes(x=BSnIPRE.gamma.bter, y=gamma)) +
	geom_point(shape=19,      # Use solid circles
			   alpha=1/2) + 	  # Make transparent
	coord_cartesian(ylim = c(-2, 11), xlim = c(-2, 6)) + 
	theme_few() +	
	theme(
			axis.text.y=element_text(size=11),
			axis.text.x=element_text(size=11)		
			
		) +
	labs(x = expression(italic("Bombus ")~gamma), y = expression(italic("Apis ")~gamma)) +
	geom_abline( slope =1, lty=2, col ="red")
p
dev.off()


#Apis vs Bombus Ortholog Gamma (Based on OrthoDB)--------------------------
pdf(file="ApisvsBombusOrthoDBGamma.pdf")
p<-ggplot(ortho, aes(x=b.gamma, y=gamma)) +
	geom_point(shape=19,      # Use solid circles
			   alpha=1/2) + 	  # Make transparent
	coord_cartesian(ylim = c(-2, 11), xlim = c(-2, 6)) + 
	theme_few() +	
	theme(
			axis.text.y=element_text(size=10),
			axis.text.x=element_text(size=10)		
			
		) +
	labs(x = expression(italic("Bombus ")~gamma), y = expression(italic("Apis ")~gamma)) + 
	geom_abline( slope =1, lty=2, col ="red")
p
dev.off()








#Gamma in Apis vs Bombus QWD genes --------------------------

#Build Plot Dataframe -----------------------
amel = merge(amel, amel.gamma, by ="GB")
degs = aggregate(degs[c(2:10)], by=list(degs$XP), mean)
names(degs)[1]="XP"

sp = rep("NDEG", nrow(degs))
sp[which(degs$FW<0.05 & degs$QG>0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "F" #F
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "Q" #Q
sp[which(degs$FW<0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG<0.05)] = "W" #W
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ>0.05 & degs$WG<0.05)] = "G"#G
degs$sp2=sp
sp = as.character(degs$sp2)
sp[sp=="F"] = "Reproductive"
sp[sp=="Q"] = "Reproductive"
sp[sp=="W"] = "Non-reproductive" #non-reproductive
sp[sp=="G"] = "Non-reproductive"
degs$sp2=sp

sp = as.character(amel$Caste)
sp[sp=="Q"] = "Reproductive"
sp[sp=="W"] = "Non-reproductive"
amel$sp2 = sp

degs.plot = degs[c(1,2,11)]; degs.plot$genus = rep("Bumblebees", nrow(degs)); names(degs.plot)[1]="GB"
amel.plot = amel[c(1,3,4)]; amel.plot$genus = rep("Honey bees", nrow(amel.plot))
qwd.plot = rbind(amel.plot, degs.plot)
qwd.plot = qwd.plot[complete.cases(qwd.plot ),] 


melted = melt(qwd.plot, id.vars=c("genus","sp2"),measure.vars=c("gamma"))
means.sem = ddply(melted, c("genus","sp2", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)



#Plot Apis and Bombus repro versus non-repro --------------------------------- 
pdf(file="QWD_rep_Selection.pdf")

p<-ggplot(means.sem, aes(x = factor(sp2), y = mean,  fill = factor(sp2))) +  
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, .9)) + 
	facet_grid(. ~ genus) +
	geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 	+
	theme_few() + 
	scale_fill_brewer(palette = "Set1") +
	theme(
		strip.text.x = element_text(size = 12.5),
		axis.title.y = element_text(size = 12.5),
		axis.text.x=element_blank(),
		axis.text=element_text(size=10),
		axis.ticks = element_blank(),
		legend.position= c(0.2,0.9)
		) +
	guides(fill=guide_legend(title=NULL)) +
	labs(
	x= "",
	y = expression(gamma)) 
p

dev.off()






# Orthology and TRG Status vs Gamma ------------------------------

#Build Plot Dataframe -----------------------
melted = melt(ortho.gamma, id.vars=c("genus","trg","n"),measure.vars=c("gamma"))
means.sem = ddply(melted, c("genus","trg", "n","variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
df = means.sem
means.sem = df[df$n=="2",] #2 supports 
trg = as.character(means.sem$trg)
trg [trg=="APOIDEA"] = "Apoidea"
trg [trg=="HYMENOPTERA"] = "Hymenoptera"
trg [trg=="INSECTA"] = "Insecta"
means.sem$trg = trg


#Plot Apis and Bombus TRG Status --------------------------------- 
pdf(file="TRGOrthologs.pdf")

p<-ggplot(means.sem, aes(x = factor(trg), y = mean,  fill = factor(trg))) + 
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, 1.5)) + 
	facet_grid(. ~ genus) +
	geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 	+
	theme_few() + 
	scale_fill_brewer(palette = "Set1") +
	theme(
		strip.text.x = element_text(size = 12.5),
		axis.title.y = element_text(size = 12.5),
		axis.text.x=element_blank(),
		axis.text=element_text(size=10),
		axis.ticks = element_blank(),
		legend.position= c(0.2,0.90)
		) +
	guides(fill=guide_legend(title=NULL)) +
	labs(
	x= "",
	y = expression(gamma)) 
	
p

dev.off()
















































































# Plot gamma relationship and histograms on side --------------------------------

	#http://www.r-bloggers.com/example-10-3-enhanced-scatterplot-with-marginal-histograms/
	#http://www.r-bloggers.com/example-8-41-scatterplot-with-marginal-histograms/

#Scatterplot/Histogram
#Define Variables and Labels
pdf(file="BombusGamma.pdf")
x=test$BSnIPRE.gamma.bter
y=test$BSnIPRE.gamma.bmel
xlab = expression(gamma~italic(" B. terrestis")~ "vs"~italic("B. impatiens"))
ylab = expression(gamma~italic(" B. melanopygus")~ "vs"~italic("B. impatiens"))
#Create Histograms and defince ranges
xhist <- hist(x, breaks=20, plot=FALSE)
yhist <- hist(y, breaks=20, plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- range(x)
yrange <- range(y)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
#Plot Center Scatterplot
par(mar=c(5,5,1,1))
plot(x, y, 
	xlim=xrange, 
	ylim=yrange, 
	xlab="", 
	ylab="",
	pch=19, 
	cex.axis=1.2,
	col=adjustcolor( "black", alpha.f = 0.1)
	)	
abline(h=mean(y),lty=2,col="grey",cex=1.2)
abline(v=mean(x),lty=2,col="grey",cex=1.2)
#Label Axes
mtext(xlab, side=1, line=2.5,adj=1, at=4,cex=1.3)
mtext(ylab, side=2, line=2.5,adj=1, at=4,cex=1.3)
#Plot Top Histogram
par(mar=c(0,3,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
#Plot Right Scatterplot
par(mar=c(3,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
dev.off()
par(def.par)

#Histogram, Bombus -----------------------

pdf(file="BombusGammaHist.pdf")
options(scipen=5)
hist(x,
	breaks=80,
	main ="",
	col = terrain.colors(70),
	xlab =  expression(gamma~" (4Nes)"),
	lwd = 1.5)
dev.off()




#Histogram, Bombus and Apis -----------------------
library(ggplot2)
test = test[c(1,21)]; names(test) = c("name","gamma")
test$genus = rep("Bumblebees", nrow(test))
names(amel.gamma) = c("name","gamma")
amel.gamma$genus = rep("Honey bees", nrow(amel.gamma))
gamma.hist = rbind(amel.gamma, test)



pdf(file="ApisBombusGammaHistogram.pdf",width=16,height=10)
ggplot(gamma.hist, aes(x=gamma)) +
geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.25,
                   colour="black", fill="white") +
geom_density(alpha=.2, fill="#FF6666") +
facet_wrap(~ genus) +
labs(
	y= "Density",
	x = expression(gamma)) +
theme_bw() + 
theme(
	strip.text.x = element_text(size = 16),
	axis.title.x = element_text(size = 16)
	)
dev.off()










# Gamma versus Alpha -----------------------------------
data <- read.table("/media/data1/bombus/Snipre/snipre_run/snipre_input_BmelToBimp.txt", header = TRUE)
#calculated alpha
test$zers = apply(test[c(3:6)],1,function(x) length(x[x=="0"]))
plot(test$BSnIPRE.gamma[test$zers<1], test$alpha[test$zers<1])

cor.test(test$BSnIPRE.gamma[test$zers<1], test$alpha[test$zers<1])
#
#        Pearson's product-moment correlation
#
#data:  test$BSnIPRE.gamma[test$zers < 1] and test$alpha[test$zers < 1]
#t = 87.479, df = 4311, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7887671 0.8102887
#sample estimates:
#      cor
#0.7997848
#



#Make a plot
pdf("AlphavsGamma.pdf")
plot(test$BSnIPRE.gamma[test$zers<1], test$alpha[test$zers<1],
	col=adjustcolor("black", alpha.f = 0.5),
	pch = 19,
	xlab = "",
	ylab = "",
	cex.axis = 1.5,
	frame.plot = F
	)
mtext(expression(gamma),side=1,line=3,cex=2.5)
mtext(expression(alpha),side=2,line=2.5,cex=2.5)
dev.off()




# QWDGF --------------------------------
#using Hollis' dataset, I'm going to look at relationships between biased genes.
	#I did this previously in BombusDEGs.r, but I'll add it here and repeat
	#In the cases where multiple gene IDs hit multiple probes, I'm going to average P values. 
	
degs = read.table(file="/media/data1/bombus/DEG/rawDEGs.txt",header=T)
degs = aggregate(degs[c(2:10)], by=list(degs$XP), mean)
names(degs)[1]="XP"
sp = rep("NDEG", nrow(degs))
sp[which(degs$FW<0.05 & degs$QG>0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "F" #F
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG>0.05)] = "Q" #Q
sp[which(degs$FW<0.05 & degs$QG<0.05 & degs$FQ<0.05 & degs$WG<0.05)] = "W" #W
sp[which(degs$FW>0.05 & degs$QG<0.05 & degs$FQ>0.05 & degs$WG<0.05)] = "G"#G
table(sp)
#   F    G NDEG    Q    W
  #10 1036 4259  216  122
degs$sp = sp
sp[sp=="NDEG"] = 1
sp[sp=="Q"] = 2
sp[sp=="G"] = 3
sp[sp=="F"] = 4
sp[sp=="W"] = 5
degs$numsp = sp
#
#sp[sp==2] = "QG"
#sp[sp==3] = "QG"
#sp[sp==4] = "WF"
#sp[sp==5] = "WF"
#degs$brsp = sp
#
#sp[sp==2] = "QF"
#sp[sp==3] = "WG"
#sp[sp==4] = "QF"
#sp[sp==5] = "WG"
#degs$repsp = sp
#

anov = aov(degs$gamma ~ degs$sp)
summary(anov)
TukeyHSD(anov)
#
#              Df Sum Sq Mean Sq F value   Pr(>F)
#degs$sp        4    8.5  2.1262   5.694 0.000143 ***
#Residuals   5638 2105.4  0.3734
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#> TukeyHSD(anov)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = degs$gamma ~ degs$sp)
#
#$`degs$sp`
#               diff           lwr         upr     p adj
#G-F    -0.005015282 -0.5348455352 0.524814972 0.9999999
#NDEG-F -0.053078433 -0.5809886240 0.474831758 0.9987602
#Q-F     0.120292077 -0.4190671746 0.659651328 0.9737891
#W-F     0.053727727 -0.4947485382 0.602203992 0.9988810
#NDEG-G -0.048063151 -0.1058261955 0.009699893 0.1546588
#Q-G     0.125307359  0.0005844589 0.250030258 0.0482717 **
#W-G     0.058743009 -0.1008615437 0.218347561 0.8535106
#Q-NDEG  0.173370510  0.0570740185 0.289667001 0.0004619 **
#W-NDEG  0.106806160 -0.0463038833 0.259916204 0.3154914
#W-Q    -0.066564350 -0.2554079831 0.122279284 0.8722218
#


#I noticed that "reproductive" genes tended to be higher, so I pulled in Hollis' data:
	#Importing Hollis' "MainRep" and "MainBrood" categories.
hol.model = read.table(file="/media/data1/bombus/DEG/hollismodel.txt",header=T)
#choose randomly
test =merge(degs, hol.model, by="XP")

anov = aov(test$gamma ~ as.factor(test$MainRep))
summary(anov)
TukeyHSD(anov)

#                           Df Sum Sq Mean Sq F value   Pr(>F)
#as.factor(test$MainRep)     2     16   7.808   19.59 3.19e-09 ***
#Residuals               13770   5487   0.399
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#> TukeyHSD(anov)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = test$gamma ~ as.factor(test$MainRep))
#
#$`as.factor(test$MainRep)`
#            diff         lwr         upr     p adj
#0--1 -0.02235219 -0.05425987 0.009555487 0.2280614
#1--1  0.05599322  0.02132285 0.090663579 0.0004523
#1-0   0.07834541  0.04883761 0.107853200 0.0000000



anov = aov(test$gamma ~ as.factor(test$MainBrood))
summary(anov)
TukeyHSD(anov)
#
#
#summary(anov)
#                             Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(test$MainBrood)     2      1  0.4609   1.154  0.316
#Residuals                 13770   5502  0.3996
#> TukeyHSD(anov)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = test$gamma ~ as.factor(test$MainBrood))
#
#$`as.factor(test$MainBrood)`
#             diff         lwr        upr     p adj
#0--1 -0.018376130 -0.04757797 0.01082571 0.3029997
#1--1 -0.006359376 -0.04158783 0.02886908 0.9060419
#1-0   0.012016754 -0.02062595 0.04465946 0.6637525
#



#Stronge effect of gamma with reproductive status but not with brood.





#Same, for AMEL ---------------------------
summary(aov(amel.plot$gamma~amel.plot$sp2))
                Df Sum Sq Mean Sq F value  Pr(>F)
amel.plot$sp2    2     24  11.976    15.8 1.6e-07 ***
Residuals     1688   1280   0.758
---
Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
10612 observations deleted due to missingness
> TukeyHSD(aov(amel.plot$gamma~amel.plot$sp2))
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = amel.plot$gamma ~ amel.plot$sp2)

$`amel.plot$sp2`
                diff         lwr          upr     p adj
NREP-NDEG  0.5002127  0.28311644  0.717309031 0.0000002
REP-NDEG   0.1839014 -0.05180573  0.419608545 0.1600526
REP-NREP  -0.3163113 -0.62806766 -0.004554994 0.0458287


# For BIMP -----------------------------------
 TukeyHSD(aov(degs.plot$gamma~degs.plot$sp2))
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = degs.plot$gamma ~ degs.plot$sp2)

$`degs.plot$sp2`
                diff         lwr       upr     p adj
NREP-NDEG 0.05425197 0.006776838 0.1017271 0.0202483
REP-NDEG  0.16804785 0.070263805 0.2658319 0.0001678
REP-NREP  0.11379589 0.009623078 0.2179687 0.0282244









#boxplot reproductive versus non-reproductive for Bombus:

amel = read.table("AMELDEGs.txt", header=T) #load in sd5 from PNAS
amel.gamma = read.table(file="/home/brock/gamma",header=T)
amel = merge(amel, amel.gamma, by ="GB")

sp = as.character(degs$sp)
sp[sp=="F"] = "REP"
sp[sp=="Q"] = "REP"
sp[sp=="W"] = "NREP"
sp[sp=="G"] = "NREP"
degs$sp2 = sp

sp = as.character(amel$Caste)
sp[sp=="Q"] = "REP"
sp[sp=="W"] = "NREP"
amel$sp2 = sp

degs.plot = degs[c(1,2,13)]; degs.plot$genus = rep("Bumblebee", nrow(degs)); names(degs.plot)[1]="GB"
amel.plot = amel[c(1,3,4)]; amel.plot$genus = rep("Honey bee", nrow(amel.plot))
qwd.plot = rbind(amel.plot, degs.plot)
qwd.plot = qwd.plot[complete.cases(qwd.plot ),]
#plot the above in 2 panel boxplots:
library(ggplot2)
 

pdf(file="QWD_rep_Selection.pdf")
ggplot(qwd.plot, aes(x = sp2, y = gamma, fill = sp2)) +
geom_boxplot(notch=TRUE, outlier.shape = 19, color="black",lwd=1.01 ) +
facet_wrap(~ genus) +
labs(x = "", y=expression(gamma))+ 
theme_bw() +
theme(
	axis.text=element_text(size=12),
	axis.ticks.x = element_blank(),
	axis.title=element_text(size=16,face="bold"),
	strip.text.x = element_text(size = 16),
	axis.text.x = element_blank(),	
	panel.border = element_rect(size=1),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
	legend.position=c(0.22,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=15),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()
	) +	
scale_fill_discrete(name="",labels=c("Non-DEG", "Non-Reproductive", "Reproductive"))

dev.off()	



# Plot gamma among QWDGF (Supplemental)--------------------------------
pdf(file="Supplemental_QWFGGamma.pdf")
boxplot(degs$gamma~degs$numsp,
	#notch = T,
	col = c("grey","gold","lightgreen","gold","lightgreen"),
		#colour is reproductive status versus not (F and Q are reproductive)
	names = c("NDEG", "Queen", "Gyne", "Foundress", "Worker"),
	ylab="",
	ylim=c(-2, 5.5),
	#names=c("Control", "Selected"),
	boxwex = 0.75,
	staplewex = 0.5,
	whisklwd = 2,
	medlwd = 2,
	whisklty = "solid",
	staplewd = 2,
	outlwd = 1,
	boxlwd = 2,
	#show.names=FALSE,
	frame=FALSE,
	cex.lab=1.2, cex=1.2, cex.axis=1.2)
mtext(expression(gamma~" (4Nes)"),side=2,line=2.2,cex=1.5)
mtext(text = c("N=216", "N=1036", "N=10", "N=122"), side = 1, line=c(2,2,2,2), at =c(2,3,4,5))
#legend("topright",inset=.0,title=NULL,c("Queen associated","Worker associated"),
#    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate(degs$gamma~degs$numsp,FUN="mean")
points(c(1,2,3,4,5),
	   c(means[1,2],means[2,2],means[3,2],means[4,2],means[5,2]),
	   pch=23,cex=2,lwd=2,bg="white")
text(c(2,3,4,5),
	 c(means[2,2]-0.07,means[3,2]-0.04,means[4,2]-0.1,means[5,2]-0.04),
     labels=format(c(means[2,2],means[3,2],means[4,2],means[5,2]),
	 format="f",
	 digits=3),
     pos=1,cex=1.2,col="black")	
segments(x0=1, x1=2, y0=5, y1=5, lty=2,lwd=2)
text(1.5, 5.2, "P = 0.00046", cex=1.2)
segments(x0=2, x1=3, y0=4.5, y1=4.5, lty=2,lwd=2)
text(2.5, 4.7, "P = 0.048", cex=1.2)
#Q-G     0.125307359  0.0005844589 0.250030258 0.0482717 **
#Q-NDEG  0.173370510  0.0570740185 0.289667001 0.0004619 **

legend("topright",inset=.0,title=NULL,c("Reproductive","Non-Reproductive"),
    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2, bty='n')
text(4.4, 4.7, "P = 0.00045", cex=1.2)

dev.off()





# Plot Woodard Main Effect Boxplots (Supplemental) --------------------------------
library(ggplot2)


pdf(file="Supplemental_ReproductiveGenesBarPlot.pdf")
genus = as.character(test$MainRep)
genus[genus=="-1"] = "A" #"Less-expressed"
genus[genus=="0"] = "B" #"No Difference"
genus[genus=="1"] = "C" #"More-expressed"
test$MainRep2 = genus
labs <- c("Less-expressed", "No Difference", "More-expressed")
abl = median(test$gamma[which(test$MainRep=="0")],na.rm=T)


ggplot(test, aes(x = as.character(MainRep2), y = gamma)) +
geom_boxplot(notch=TRUE, outlier.shape = 19, color="black",lwd=1.01 ) +
labs(x = "", y=expression(gamma))+ 
theme_bw() +
scale_x_discrete(labels= labs) +  
geom_abline(aes(slope=0,intercept=abl,color="red")) +
theme(
	axis.text.x  = element_text(vjust=0.5, size=14)
	)

dev.off()




pdf(file="Supplemental_BroodCareGenesBarPlot.pdf")
genus = as.character(test$MainBrood)
genus[genus=="-1"] = "A" #"Less-expressed"
genus[genus=="0"] = "B" #"No Difference"
genus[genus=="1"] = "C" #"More-expressed"
test$MainBrood2 = genus
labs <- c("Less-expressed", "No Difference", "More-expressed")
abl = median(test$gamma[which(test$MainBrood=="0")],na.rm=T)


ggplot(test, aes(x = as.character(MainBrood2), y = gamma)) +
geom_boxplot(notch=TRUE, outlier.shape = 19, color="black",lwd=1.01 ) +
labs(x = "", y=expression(gamma))+ 
theme_bw() +
scale_x_discrete(labels= labs) +  
geom_abline(aes(slope=0,intercept=abl,color="red")) +
theme(
	axis.text.x  = element_text(vjust=0.5, size=14)
	)

dev.off()




























#AZ wants the same plot above but as a bar plot. OK. 

library(ggplot2)
library(plyr)
library(reshape2)

melted = melt(degs, id.vars="numsp",measure.vars=c("gamma"))
means.sem = ddply(melted, c("numsp", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)

pdf(file="AZQGFWBarPlot.pdf")
ggplot(means.sem, aes(x = factor(numsp), y = mean)) +  
  geom_bar(position = "dodge",stat="identity",  fill="lightblue",color="black" ) + 
   theme(axis.title.x = element_blank(), axis.text.x  = element_text(angle=0, vjust=0.5, size=16)) +  
  ylab("Gamma") +
 scale_x_discrete(labels=c("NDEG", "Queen", "Gyne", "Foundress", "Worker")) +
 geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 

 dev.off()


#AZ also wants a barplot by "reproductive" and "non-reproductive".
#used "hol.model" merge above.

melted = melt(test, id.vars="MainRep",measure.vars=c("gamma"))
means.sem = ddply(melted, c("MainRep", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)


pdf(file="ReproductiveGenesBarPlot.pdf")
ggplot(means.sem, aes(x = factor(MainRep), y = mean)) +  
  geom_bar(position = "dodge",stat="identity",  fill="lightblue",color="black" ) + 
   theme(axis.title.x = element_blank(), axis.text.x  = element_text(angle=0, vjust=0.5, size=16)) +  
  ylab("Gamma") +
scale_x_discrete(labels=c("Under-expressed", "NDEG", "Over-expressed")) +
 geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 
 dev.off()


 
melted = melt(test, id.vars="MainBrood",measure.vars=c("gamma"))
means.sem = ddply(melted, c("MainBrood", "variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
 
  
pdf(file="BroodCareGenesBarPlot.pdf")
ggplot(means.sem, aes(x = factor(MainBrood), y = mean)) +  
  geom_bar(position = "dodge",stat="identity",  fill="lightblue",color="black" ) + 
   theme(axis.title.x = element_blank(), axis.text.x  = element_text(angle=0, vjust=0.5, size=16)) +  
  ylab("Gamma") +
scale_x_discrete(labels=c("Under-expressed", "NDEG", "Over-expressed")) +
 geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 
 dev.off()
 
 
 
 
 


test =merge(degs, hol.model, by="XP")

anov = aov(test$gamma ~ as.factor(test$MainRep))
summary(anov)
TukeyHSD(anov)

#                           Df Sum Sq Mean Sq F value   Pr(>F)
#as.factor(test$MainRep)     2     16   7.808   19.59 3.19e-09 ***
#Residuals               13770   5487   0.399
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#> TukeyHSD(anov)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = test$gamma ~ as.factor(test$MainRep))
#
#$`as.factor(test$MainRep)`
#            diff         lwr         upr     p adj
#0--1 -0.02235219 -0.05425987 0.009555487 0.2280614
#1--1  0.05599322  0.02132285 0.090663579 0.0004523
#1-0   0.07834541  0.04883761 0.107853200 0.0000000








# GO Analyses -----------------------------------------------------------------
	#GO terms for highest gamma genes?
	
high.genes = as.character(test$Group[test$BSnIPRE.gamma.bter > 1])
bimptoFBGN = read.table(file="/media/data1/bombus/blast/BestDMELtoBIMPout.txt",header=F)
write.list(bimptoFBGN$V2, file="BIMPUniverse")
write.list(bimptoFBGN$V2[bimptoFBGN$V1 %in% high.genes], file="BIMPGenes")






# Prepare gene universe and your test set --------------------------------
#GOstatsBEE.r
# Libraries and functions ------------------
source("http://bioconductor.org/biocLite.R")
#biocLite("GOstats")
library("GOstats")
library("AnnotationForge")
library("GSEABase")
#available.dbschemas() #lists available data packages
library("drosophila2.db") #Download a database from those that are available
#library("KEGG.db")


#universe
	#= goframeData$frame.gene_id #define your gene universe
universe = read.table(file="BIMPUniverse",header=F)
universe = unlist(as.character(universe$V1[!is.na(universe$V1)]))
universe = universe[!duplicated(universe)]
#genes
genes = as.character(unique(unlist(read.table(file="BIMPGenes")$V1)))




#load database for GO analyses and prepare GO data.frame and GeneSetCollection ------
k = (keys(drosophila2.db,keytype="GO")) # I read in our dros. orthologs
orths = universe
goframeData = select(drosophila2.db, keys=k, cols=c("EVIDENCE", "FLYBASE"), keytype="GO")
goframeData = goframeData[goframeData$FLYBASE %in% orths, ]
goframeData = goframeData[!duplicated(goframeData$FLYBASE), ]
names(goframeData) =c("frame.go_id", "frame.Evidence","frame.ontology", "frame.gene_id")
goframeData$frame.ontology = NULL
#goframeData.BP = goframeData[goframeData$frame.ontology =="BP",]
#goframeData.BP$frame.ontology = NULL
goFrame = GOFrame(goframeData,organism="Bombus impatiens")
goAllFrame = GOAllFrame(goFrame) #this is dropping my geneIDs.
gsc = GeneSetCollection(goAllFrame, setType = GOCollection())
	


#  Perform tests ----------------
#Kindly provded by Sasha (thanks, man)
outcomeBP <- hyperGTest(GSEAGOHyperGParams(name = "lrt", 
			geneSetCollection=gsc,
			geneIds = as.character(genes),
			universeGeneIds = as.character(universe),
			ontology = "BP",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over"))
outcomeMF <- hyperGTest(GSEAGOHyperGParams(name = "lrt", 
			geneSetCollection=gsc,
			geneIds = as.character(genes),
			universeGeneIds = as.character(universe),
			ontology = "MF",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over"))
outcomeCC <- hyperGTest(GSEAGOHyperGParams(name = "lrt", 
			geneSetCollection=gsc,
			geneIds = as.character(genes),
			universeGeneIds = as.character(universe),
			ontology = "CC",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over"))
#head(summary(outcomeBP))
#head(summary(outcomeMF))			
#head(summary(outcomeCC))
		








# APIS versus BOMBUS selection -----------------------------------------------------------------

#With BLAST-Best-Matches
conv = read.table(file="/media/data1/bombus/GFF/BIMP_XP_GBoldnew_FBGN_recip.txt",header=F)
test = merge(gm.bter, gm.bmel, by="Group",suffix=c(".bter",".bmel")) 
test = test[c(1,21,10)]
names(test)[1] = names(conv)[1] = "XP"
test = merge(test, conv, by = "XP")
test = test[complete.cases(test),]
names(test)[5]="GB"
amel.gamma = read.table(file="/home/brock/gamma",header=T)
test = merge(test, amel.gamma, by = "GB")


 
cor.test(test$BSnIPRE.gamma.bter, test$gamma) 
#
#data:  test$BSnIPRE.gamma.bter and test$gamma
#t = 15.118, df = 7165, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1532929 0.1981667
#sample estimates:
#      cor
#0.1758212

cor.test(test$BSnIPRE.gamma.bter[test$gamma>1 & test$BSnIPRE.gamma.bter>1], test$gamma[test$gamma>1 & test$BSnIPRE.gamma.bter>1])
#
#data:  test$BSnIPRE.gamma.bter[test$gamma > 1 & test$BSnIPRE.gamma.bter >  and test$gamma[test$gamma > 1 & test$BSnIPRE.gamma.bter > 1]    1] and test$gamma[test$gamma > 1 & test$BSnIPRE.gamma.bter > 1]
#t = -0.74124, df = 218, p-value = 0.4593
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.18120947  0.08267994
#sample estimates:
#        cor
#-0.05013983
#


cor.test(test$BSnIPRE.gamma.bter[test$gamma>1 | test$BSnIPRE.gamma.bter>1], test$gamma[test$gamma>1 | test$BSnIPRE.gamma.bter>1])
#
#data:  test$BSnIPRE.gamma.bter[test$gamma > 1 | test$BSnIPRE.gamma.bter >  and test$gamma[test$gamma > 1 | test$BSnIPRE.gamma.bter > 1]    1] and test$gamma[test$gamma > 1 | test$BSnIPRE.gamma.bter > 1]
#t = -14.523, df = 1658, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.3779354 -0.2925484
#sample estimates:
#      cor
#-0.335932
#






ap.out = quantile(amel.gamma$gamma,0.95)
bomb.out = quantile(test$BSnIPRE.gamma.bter,0.95)
nrow(test[test$gamma>ap.out & test$BSnIPRE.gamma.bter>bomb.out ,])

write.list(amel.gamma$GB[amel.gamma$gamma > ap.out],file="highamel")



















# Plot Gamma in Bombus and Alpha for BBM ------------------

pdf(file="ApisvsBombusBBMGamma.pdf")
ggplot(test, aes(x=BSnIPRE.gamma.bter, y=gamma)) +
    geom_point(shape=19,      # Use solid circles
               alpha=1/5) + 	  # Make transparent
labs(x = expression(italic("Bombus ")~gamma), y = expression(italic("Apis ")~gamma)) +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	#axis.ticks.x = element_blank(),
	axis.title=element_text(size=16,face="bold"),
	#axis.text.x = element_blank(),	
	panel.border = element_rect(size=1)
	) 

dev.off()



# Corr.Plot Gamma in Bombus and Apis for BBM ------------------
df.plot = c()
for(i in seq(-1,2.4,.2)){
	x = cor.test(test$BSnIPRE.gamma.bter[test$BSnIPRE.gamma.bter > i & test$gamma  > i], test$gamma[test$BSnIPRE.gamma.bter > i & test$gamma  > i])
	print(c(i, x$p.value, as.numeric(x$estimate)))
	df.plot = rbind(df.plot,c(i, x$p.value, as.numeric(x$estimate)) )
}

df.plot = data.frame(df.plot)
cols = as.numeric(df.plot$X2)
cols[cols<0.05] = "P < 0.05" ;  cols[cols!="P < 0.05"] = "P > 0.05"


pdf(file="ApisGammaCorrelationPlot.pdf")
ggplot(df.plot, aes(x=X1, y=X3, col = cols)) +
geom_point(shape=19) + 	  
labs(x = expression(gamma ~"Cutoff") , y = "Pearson Correlation Coefficient (r)") +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	axis.title=element_text(size=14,face="bold"),
	panel.border = element_rect(size=1),
	legend.position=c(0.20,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=13),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()	
	) 
dev.off()
	
	
	
	
	
# Corr.Plot Gamma in Bombus ONLY  for BBM ------------------
df.plot = c()
for(i in seq(-1,2.4,.2)){
	x = cor.test(test$BSnIPRE.gamma.bter[test$BSnIPRE.gamma.bter > i ], test$gamma[test$BSnIPRE.gamma.bter > i ])
	print(c(i, x$p.value, as.numeric(x$estimate)))
	df.plot = rbind(df.plot,c(i, x$p.value, as.numeric(x$estimate)) )
}

df.plot = data.frame(df.plot)
cols = as.numeric(df.plot$X2)
cols[cols<0.05] = "P < 0.05" ;  cols[cols!="P < 0.05"] = "P > 0.05"


pdf(file="BombusONLYGammaCorrelationPlot.pdf")
ggplot(df.plot, aes(x=X1, y=X3, col = cols)) +
geom_point(shape=19) + 	  
labs(x = expression(gamma ~"Cutoff") , y = "Pearson Correlation Coefficient (r)") +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	axis.title=element_text(size=14,face="bold"),
	panel.border = element_rect(size=1),
	legend.position=c(0.20,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=13),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()	
	) 
dev.off()

	
	
	
# Proportion of highGamma genes in Bombus ONLY  for BBM for BINS ------------------
df.plot = c()
bin.size = 0.5
for(i in seq(0.5,5,bin.size)){
	#window
	high.genes = length(test$gamma[which(test$BSnIPRE.gamma.bter > i - bin.size & test$BSnIPRE.gamma.bter < i + bin.size & test$gamma > 1)])
	all.genes = length(test$gamma[which(test$BSnIPRE.gamma.bter > i - bin.size & test$BSnIPRE.gamma.bter < i + bin.size )])
	#all.else
	all.high.genes = length(test$gamma[which(test$BSnIPRE.gamma.bter < i - bin.size | test$BSnIPRE.gamma.bter > i+ bin.size & test$gamma > 1)])
	all.all.genes = length(test$gamma[which(test$BSnIPRE.gamma.bter < i - bin.size | test$BSnIPRE.gamma.bter > i+ bin.size)])
	matr = matrix(nc=2,nr=2)
	matr[1,] =  c(high.genes,all.genes )
	matr[2,] = c(all.high.genes,all.all.genes )
	x = chisq.test(matr, )
	prop = 	(high.genes/all.genes) / (all.high.genes/all.all.genes)
	prop2 = 	(high.genes/all.genes)
	#i = mean(c(i - bin.size, i+bin.size))
	print(c(i, x$p.value, prop))
	df.plot = rbind(df.plot,c(i, x$p.value, prop, prop2) )
}

df.plot = data.frame(df.plot)
cols = as.numeric(df.plot$X2)
cols[cols<0.05] = "P < 0.05" ;  cols[cols!="P < 0.05"] = "P > 0.05"

pdf(file="RelativeProportionOutlierBinned.pdf")
ggplot(df.plot, aes(x=X1, y=X4, col=cols)) +
geom_point(shape=19) + 	  
labs(x = expression(gamma ~"Cutoff") , y = "Relative Proportion Outlier") +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	axis.title=element_text(size=14,face="bold"),
	panel.border = element_rect(size=1),
	legend.position=c(0.20,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=13),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()	
	) 
dev.off()


#######
with(test,plot(BSnIPRE.gamma.bter,gamma))
abline(h=1, v=1, lty=2,col="red")
#########
abline(v=c(1.5, 2), lty=2,col="blue")








# Proportion of highGamma genes in Apis ONLY  for BBM for BINS ------------------
df.plot = c()
bin.size = 0.5
for(i in seq(0.5,5,bin.size)){
	#window
	high.genes = length(test$BSnIPRE.gamma.bter[which(test$gamma > i - bin.size & test$gamma < i + bin.size & test$BSnIPRE.gamma.bter > 1)])
	all.genes = length(test$BSnIPRE.gamma.bter[which(test$gamma > i - bin.size & test$gamma < i + bin.size )])
	#all.else
	all.high.genes = length(test$BSnIPRE.gamma.bter[which(test$gamma < i - bin.size | test$gamma > i+ bin.size & test$BSnIPRE.gamma.bter > 1)])
	all.all.genes = length(test$BSnIPRE.gamma.bter[which(test$gamma < i - bin.size | test$gamma > i+ bin.size)])
	matr = matrix(nc=2,nr=2)
	matr[1,] =  c(high.genes,all.genes )
	matr[2,] = c(all.high.genes,all.all.genes )
	x = chisq.test(matr, )
	prop = 	(high.genes/all.genes) / (all.high.genes/all.all.genes)
	#i = mean(c(i - bin.size, i+bin.size))
	print(c(i, x$p.value, prop))
	df.plot = rbind(df.plot,c(i, x$p.value, prop) )
}

df.plot = data.frame(df.plot)
cols = as.numeric(df.plot$X2)
cols[cols<0.05] = "P < 0.05" ;  cols[cols!="P < 0.05"] = "P > 0.05"

pdf(file="RelativeProportionOutlierBinned.pdf")
ggplot(df.plot, aes(x=X1, y=X3, col=cols)) +
geom_point(shape=19) + 	  
labs(x = expression(gamma ~"Cutoff") , y = "Relative Proportion Outlier") +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	axis.title=element_text(size=14,face="bold"),
	panel.border = element_rect(size=1),
	legend.position=c(0.20,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=13),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()	
	) 
dev.off()


#######
with(test,plot(BSnIPRE.gamma.bter,gamma))
abline(h=1, v=1, lty=2,col="red")
#########
abline(v=c(1.5, 2), lty=2,col="blue")







# Proportion of highGamma genes in Bombus ONLY  for ORTHO for BINS ------------------



df.plot = c()
bin.size = 0.5
for(i in seq(0.5,5,bin.size)){
	#window
	high.genes = length(ortho$ortho.gamma[which(ortho$b.gamma > i - bin.size & ortho$b.gamma < i + bin.size & ortho$ortho.gamma > 1)])
	all.genes = length(ortho$ortho.gamma[which(ortho$b.gamma > i - bin.size & ortho$b.gamma < i + bin.size )])
	#all.else
	all.high.genes = length(ortho$ortho.gamma[which(ortho$b.gamma < i - bin.size | ortho$b.gamma > i+ bin.size & ortho$ortho.gamma > 1)])
	all.all.genes = length(ortho$ortho.gamma[which(ortho$b.gamma < i - bin.size | ortho$b.gamma > i+ bin.size)])
	matr = matrix(nc=2,nr=2)
	matr[1,] =  c(high.genes,all.genes )
	matr[2,] = c(all.high.genes,all.all.genes )
	x = chisq.test(matr )
	prop = 	(high.genes/all.genes) / (all.high.genes/all.all.genes)
	prop2 = 	(high.genes/all.genes)
	#i = mean(c(i - bin.size, i+bin.size))
	print(c(i, x$p.value, prop))
	df.plot = rbind(df.plot,c(i, x$p.value, prop, prop2) )
}

df.plot = data.frame(df.plot)
cols = as.numeric(df.plot$X2)
cols[cols<0.05] = "P < 0.05" ;  cols[cols!="P < 0.05"] = "P > 0.05"

pdf(file="RelativeProportionOutlierBinnedORTHO.pdf")
ggplot(df.plot, aes(x=X1, y=X4, col=cols)) +
geom_point(shape=19) + 	  
labs(x = expression(gamma ~"Cutoff") , y = "Relative Proportion Outlier") +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	axis.title=element_text(size=14,face="bold"),
	panel.border = element_rect(size=1),
	legend.position=c(0.20,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=13),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()	
	) 
dev.off()


#######
with(ortho,plot(b.gamma,ortho.gamma))
abline(h=1, v=1, lty=2,col="red")
#########
abline(v=c(1.5, 2), lty=2,col="blue")






































# Plot Gamma in Bombus and Alpha for OrthoDB ------------------
	
pdf(file="ApisvsBombusOrthoDBGamma.pdf")
ggplot(ortho, aes(x=b.gamma, y=ortho.gamma)) +
    geom_point(shape=19,      # Use solid circles
               alpha=1/5) + 	  # Make transparent
labs(x = expression(italic("Bombus ")~gamma), y = expression(italic("Apis ")~gamma)) +
theme_bw() +	
theme(
	axis.text=element_text(size=12),
	#axis.ticks.x = element_blank(),
	axis.title=element_text(size=16,face="bold"),
	#axis.text.x = element_blank(),	
	panel.border = element_rect(size=1)
	) 

dev.off()








#With OrthoDB
for(i in seq(-1,5,.5)){
	x = cor.test(ortho$b.gamma[ortho$b.gamma > i &ortho$ortho.gamma  > i],  ortho$ortho.gamma[ortho$b.gamma > i & ortho$ortho.gamma  > i])
	print(c(i, x$p.value, as.numeric(x$estimate)))

}







 
 
 
 
#GO Term overlap?
	#I ran DAVID with ONLY GO terms for gamma>1 in AMEL and BIMP. (BIMPGO and AMELGO).

GO = read.table(file="/media/data1/bombus/workingFiles/GOOverlap.txt",header=T)
write.table(table(GO), file="/media/data1/bombus/workingFiles/GOOverlapTABLE.txt")
GO = read.table(file="/media/data1/bombus/workingFiles/GOOverlapTABLE.txt")


table(aggregate(GO$GO, by=list(GO$GO), length)$x)
#Overlap is very very slight. Same networks not acted on by selection. 
#  1   2
#256   4
goover = (aggregate(GO$GO, by=list(GO$GO), length))
goover[goover$x=="2",]

#       Group.1 x
#GO:0016021	integral to membrane
#GO:0031224	intrinsic to membrane
#GO:0050877	neurological system process
#GO:0051606	detection of stimulus



# Venn Diagram --------------------------
library(VennDiagram)

ca = nrow((GO[rowSums(GO)=="2",])) #cross-area
a1 = nrow((GO[rowSums(GO[1])=="1",])) - ca #Bombus-specific
a2 = nrow((GO[rowSums(GO[2])=="1",])) - ca#Apis Specific

grid.newpage()
draw.pairwise.venn(area1 = a1, area2 = a2, cross.area = ca, category = c("Dog People", 
    "Cat People"))




























# Orthology and TRG Status ------------------------------

#load ortho list and bombus gamma.
ortho = read.table(file="/media/data1/bombus/workingFiles/OrthoGammaBIMPvAPIS.txt",header=T)
#headers: 1 - GB or BIMP v1, class , orthoDB OGID, the orthologous GB or BIMP ID, gamma for bombus, normalized bombus gamma, apis gamma, normalized apis gamma, gamma of the ortholog, normalized ortholog gamma

ortho.high = ortho[which(ortho$ortho.gamma >1 | ortho$b.gamma > 1),]
with(ortho.high, plot(ortho.gamma, b.gamma))
with(ortho.high, cor.test(ortho.gamma, b.gamma))

plot(ortho$b.gamma,ortho$ortho.gamma) 
plot(ortho$b.gamma[ortho$class=="INSECTA"],ortho$ortho.gamma[ortho$class=="INSECTA"], pch=19)
reg1 = lm(ortho$b.gamma[ortho$class=="INSECTA"]~ortho$ortho.gamma[ortho$class=="INSECTA"])
abline(reg1)
points(ortho$b.gamma[ortho$class=="HYMENOPTERA"],ortho$ortho.gamma[ortho$class=="HYMENOPTERA"], col="grey", pch=19) 
reg1 = lm(ortho$b.gamma[ortho$class=="HYMENOPTERA"]~ortho$ortho.gamma[ortho$class=="HYMENOPTERA"])
abline(reg1,col="grey")
points(ortho$b.gamma[ortho$class=="APOIDEA"],ortho$ortho.gamma[ortho$class=="APOIDEA"],col="red", pch=19) 
reg1 = lm(ortho$b.gamma[ortho$class=="APOIDEA"]~ortho$ortho.gamma[ortho$class=="APOIDEA"])
abline(reg1,col="red", lty=2)
#No difference, really in correlation of any of these.


summary(aov(orth2$b.gamma~orth2$class) )
#              Df Sum Sq Mean Sq F value   Pr(>F)
#orth2$class    3      7  2.2704    5.98 0.000456 ***
#Residuals   8457   3211  0.3796
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#3594 observations deleted due to missingness
TukeyHSD(aov(orth2$b.gamma~orth2$class) )
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = orth2$b.gamma ~ orth2$class)
#
#$`orth2$class`
#                           diff        lwr         upr     p adj
#BOMBUS-APOIDEA      -0.04698990 -0.3333905  0.23941071 0.9748085
#HYMENOPTERA-APOIDEA -0.14715205 -0.3067481  0.01244402 0.0831583
#INSECTA-APOIDEA     -0.19942204 -0.3450645 -0.05377952 0.0024633
#HYMENOPTERA-BOMBUS  -0.10016215 -0.3565200  0.15619571 0.7471011
#INSECTA-BOMBUS      -0.15243213 -0.4003438  0.09547952 0.3901599
#INSECTA-HYMENOPTERA -0.05226999 -0.1223113  0.01777129 0.2206330
#

TukeyHSD(aov(orth2$ortho.gamma~orth2$class) )
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = orth2$ortho.gamma ~ orth2$class)
#
#$`orth2$class`
#                          diff        lwr         upr     p adj
#HYMENOPTERA-APOIDEA -0.2593090 -0.4442894 -0.07432858 0.0029306
#INSECTA-APOIDEA     -0.5978388 -0.7691250 -0.42655266 0.0000000
#INSECTA-HYMENOPTERA -0.3385298 -0.4134953 -0.26356433 0.0000000
#
summary(aov(orth2$ortho.gamma~orth2$class) )
#              Df Sum Sq Mean Sq F value Pr(>F)
#orth2$class    2     99   49.57   87.07 <2e-16 ***
#Residuals   9145   5207    0.57
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#2907 observations deleted due to missingness
#
#







#plot the above in 2 panel boxplots:
library(ggplot2)
 
genus = as.character(ortho.plot$genus)
genus[genus=="APIS"] = "Honey bees"
genus[genus=="BOMBUS"] = "Bumblebees"
ortho.plot$genus2=genus
pdf(file="TRGOrthologs.pdf")
ggplot(ortho.plot, aes(x = class, y = b.gamma, fill = class)) +
geom_boxplot(notch=TRUE, outlier.shape = 19, color="black",lwd=1.01 ) +
facet_wrap(~ genus2) +
labs(x = "", y=expression(gamma))+ 
theme_bw() +
theme(
	axis.text=element_text(size=12),
	axis.ticks.x = element_blank(),
	axis.title=element_text(size=16,face="bold"),
	strip.text.x = element_text(size = 16),
	axis.text.x = element_blank(),	
	panel.border = element_rect(size=1),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
	legend.position=c(0.20,0.75),
	legend.title=element_blank(),
	legend.text=element_text(size=15),
	#legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
	legend.key = element_blank()
	) +
scale_fill_discrete(name="",labels=c("Apoidea", "Hymenoptera", "Insecta"))
dev.off()	




#plot Gamma vs Gamma for Bombus and Apis genes:
ortho = ortho.plot
ortho =  ortho[!duplicated(ortho$BvAortholog),]



pdf(file="BombusApisGamma.pdf")
cols = rep("black", nrow(ortho))
#cols[test$gamma>1 | test$BSnIPRE.gamma.bter>1] ="lightred"
#cols[test$gamma>1 & test$BSnIPRE.gamma.bter>1] ="red"
plot(ortho$b.gamma, ortho$ortho.gamma,
	pch=19,
	col=adjustcolor(cols, alpha.f = 0.2),
	xlab = expression(italic("Bombus ")~gamma),
	ylab = expression(italic("Apis ")~gamma)
	)


plot(ortho$b.gamma[ortho$b.gamma>1 & ortho$ortho.gamma>1], ortho$ortho.gamma[ortho$b.gamma>1 & ortho$ortho.gamma>1])
cor.test(ortho$b.gamma, ortho$ortho.gamma)
plot(ortho$b.gamma[ortho$class!="APOIDEA"], ortho$ortho.gamma[ortho$class!="APOIDEA"])



# More overlap than chance? 
	#		BIMP
	#			>1	<1
	#APIS	>1  267	448
	#		<1	989	5467
	#
	#The two-tailed P value is less than 0.0001




#Genes with high gamma in one lineage don't usually in another (BIMP versus 
test1=test[test$gamma>1 | test$BSnIPRE.gamma.bter>1,]
plot(test1$BSnIPRE.gamma.bter, test1$gamma)
x = glm(test1$BSnIPRE.gamma.bter~test1$gamma)
Y = predict(x,test1)
 




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 









# TF analysis-----------------------------------------------------------------
	#I took all TFs from FlyBase (go-id:GO:0003700) for DMEL only
test = merge(gm.bter, gm.bmel, by="Group",suffix=c(".bter",".bmel"))
test = test[c(1,21,45)]
tfs = read.table(file="FlyBase_IDs.txt", header=F); tfs$TF = rep("TF", nrow(tfs))
names(tfs)[1]="FBGN"
bimptoFBGN = read.table(file="/media/data1/bombus/blast/BestDMELtoBIMPout.txt",header=F)
names(bimptoFBGN)[1]="Group"
test = merge(test, bimptoFBGN, all.x=T,by="Group")
test = test[c(1,2,3,4)];names(test)[4]="FBGN"
test = merge(test, tfs, all.x=T,by="FBGN")
tf.status = as.character(test$TF); tf.status[is.na(tf.status)]="non"
test$TF = tf.status
summary(aov(test$BSnIPRE.gamma.bter~test$TF))
#               Df Sum Sq Mean Sq F value   Pr(>F)
#test$TF         1      5   4.747    11.6 0.000662 ***
#Residuals   10006   4095   0.409
	#lower TF versus genome.

	
summary(aov(test$BSnIPRE.gamma.bter[test$BSnIPRE.gamma.bter>1]~test$TF[test$BSnIPRE.gamma.bter>1]))	
#no difference for gamma>1	
	

#AMEL versus BIMP TF 
FBGN=read.table(file="/media/data1/forty3/brock/scripts/FBGN.txt",header=T)
names(FBGN)[2]="FBGN"
amel.test = merge(amel.gamma, FBGN,  by="GB",all.x=T)
amel.test = merge(amel.test, tfs,  by="FBGN",all.x=T)
tf.status = as.character(amel.test$TF); tf.status[is.na(tf.status)]="non"
amel.test$TF = tf.status

summary(aov(amel.test$gamma~amel.test$TF))
#                Df Sum Sq Mean Sq F value Pr(>F)
#amel.test$TF     1      0  0.0126   0.024  0.878
#Residuals    12301   6558  0.5331
#
#no difference in AMEL



summary(aov(amel.test$gamma[amel.test$gamma>1]~amel.test$TF[amel.test$gamma>1]))
#slightly different, TF lower.






	
	


##looking for differentially expressed genes between W Q G F categories in the microarray data (Woodard et all 2014)

a = read.table("/home/jen/diffex.txt", head = T)

####################################################################################
# This part removes duplicate probes and takes average p for genes expressed when
# they are targets of more than one probe

b = a[!duplicated(a$Seq),]
WF <-aggregate(b$WF, by=list(Gene=b$Gene), FUN=mean)
WF = WF[which(WF$x<0.05),]

WG <-aggregate(b$WG, by=list(Gene=b$Gene), FUN=mean)
WG = WG[which(WG$x<0.05),]

QG <-aggregate(b$QG, by=list(Gene=b$Gene), FUN=mean)
QG = QG[which(QG$x<0.05),]

QF <-aggregate(b$QF, by=list(Gene=b$Gene), FUN=mean)
QF = QF[which(QF$x<0.05),]
####################################################################################

WF = WF[,c(1,2,3)]
WG = WG[,c(1,2,6)]
QG = QG[,c(1,2,4)]
QF = QF[,c(1,2,5)]

aa = merge(WF,WG, by ="Gene")
W<-aa[which(!aa$Gene%in%QG$Gene),]
W<-W[which(!W$Gene%in%QF$Gene),]

aa = merge(WF,QF, by ="Gene")
F<-aa[which(!aa$Gene%in%WG$Gene),]
F<-F[which(!F$Gene%in%QG$Gene),]

aa = merge(QG,WG, by ="Gene")
G<-aa[which(!aa$Gene%in%WF$Gene),]
G<-G[which(!G$Gene%in%QF$Gene),]

aa = merge(QG,QF, by ="Gene")
Q<-aa[which(!aa$Gene%in%WF$Gene),]
Q<-Q[which(!Q$Gene%in%WG$Gene),]

# four data frames(Q,G,F,W), each with only significant expression values

columns <- c("Gene", "Gamma")
gm = read.table("/media/data1/bombus/Snipre/allgamma.txt", head=FALSE,skip=1,col.names=columns)

Wg = merge(W,gm, by ="Gene")
Fg = merge(F,gm, by ="Gene")
Qg = merge(Q,gm, by ="Gene")
Gg = merge(G,gm, by ="Gene")

#note we lose a few XPs above because of the Snipre filters
#we parse these four files together
Wg$Class <- "W"
Fg$Class <- "F"
Qg$Class <- "Q"
Gg$Class <- "G"

List = rbind(Wg,Fg,Qg,Gg)
#write.table(List, file="/home/jen/Bter_gamma_expression.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
List$Gene <- List$Gene.x

plot1 = List[,c(5,4)]

boxplot(Gamma~Class,data=plot1, main="Expressed Gamma", notch=T, xlab="Class", ylab="Gamma")
bb=aov(Gamma~Class,data=plot1)
summary(bb)
TukeyHSD(bb)

# if you want to take out foundresses (low n) or gynes
plot2<-plot1[which(!plot1$Class=="F"),]
plot2<-plot2[which(!plot2$Class=="G"),]

boxplot(x~Class,data=plot2, main="Expressed Gamma", notch=T, xlab="Class", ylab="Gamma")
aa=aov(x~Class,data=plot2)
summary(aa)
TukeyHSD(aa)

#### looking at differential expression in reproductive and brood caring individuals
#### this time it is from the pvalues Woodards two columns for brood and reproductive
a = read.table("/home/jen/diffex2.txt", head = T)
b = a[!duplicated(a$Seq),]

repl <-aggregate(b$Rep, by=list(Gene=b$Gene), FUN=mean)
repl = repl[which(repl$x<0.05),]

brood <-aggregate(b$Brood, by=list(Gene=b$Gene), FUN=mean)
brood = brood[which(brood$x<0.05),]

columns <- c("Gene", "Gamma")
gm = read.table("/media/data1/bombus/Snipre/allgamma.txt", head=FALSE,skip=1,col.names=columns)

Re = merge(repl,gm, by ="Gene")
B = merge(brood,gm, by ="Gene")

Re$Class<-"R"
B$Class<-"B"

List = rbind(Re,B)

plot1 = List[,c(4,3)]

boxplot(Gamma~Class,data=plot1, main="Expressed Gamma", notch=T, xlab="Class", ylab="Gamma")
bb=aov(Gamma~Class,data=plot1)
summary(bb)
TukeyHSD(bb)
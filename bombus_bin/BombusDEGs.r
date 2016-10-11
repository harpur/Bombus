

names(gbs) = c("XP", "V2", "GB")
deg = merge(deg, gbs, by="GB")


deg = deg[!duplicated(deg$XP),] # Fixing duplicate XPs 
sp = rep("NDEG", nrow(deg))
sp[which(deg$FW<0.05 & deg$QG>0.05 & deg$FQ<0.05 & deg$WG>0.05)] = "F" #F
sp[which(deg$FW>0.05 & deg$QG<0.05 & deg$FQ<0.05 & deg$WG>0.05)] = "Q" #Q
sp[which(deg$FW<0.05 & deg$QG<0.05 & deg$FQ<0.05 & deg$WG<0.05)] = "W" #W
sp[which(deg$FW>0.05 & deg$QG<0.05 & deg$FQ>0.05 & deg$WG<0.05)] = "G"#G
deg$sp = sp




boxplot(deg$gamma~deg$sp)

boxplot(deg$gamma~deg[,16])


#to repeat this, need GBold dataset, SNIPRE output, and Jen's conversion table. Gamma.


deg1 = deg[deg$sp!="F",]
deg1 = deg1[deg1$sp!="G",]
x=(aov(deg1$gamma~as.factor(deg1$sp)))

#TukeyHSD(x)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = deg1$gamma ~ as.factor(deg1$sp))
#
#$`as.factor(deg1$sp)`
#             diff          lwr       upr     p adj
#Q-NDEG 0.09499531  0.007286688 0.1827039 0.0299647
#W-NDEG 0.10694123  0.010780920 0.2031015 0.0248256
#W-Q    0.01194591 -0.114050889 0.1379427 0.9731291
#

summary(x)
                     Df Sum Sq Mean Sq F value  Pr(>F)   
#as.factor(deg1$sp)    2    4.6  2.3212    6.23 0.00199 **
#Residuals          4366 1626.7  0.3726               


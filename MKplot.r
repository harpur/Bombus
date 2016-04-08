#!/usr/bin/Rscript
###
# Plot MK results
###



#load source ---------------------------------------
source(file="/GenomeR/VarFunct.r")


# Load dataframe -----------------------------------
mk = read.table(file = "MKresults", header=T)

# Assess significance cut-offs ---------------------
p = 0.05
x = FDRcontrol(mk$MKp,0.1)
np = range(mk$MKp[x])[2]

# Plot result  ---------------------
pdf("NI.pdf")
plot(-log10(mk$NI), 
	-log10(mk$MKp),	
	xlab = "-Log10(MK Neutrality Index)",
	ylab = "-Log10(MK P-value)",
	pch = 19
	)
abline(h=-log10(p), col="red", lty = 2 )
abline(h=-log10(np), col="red", lty = 2 )
dev.off()
## @knitr joineQTLdatasets

t <- unlist(strsplit(tissue, split=""))[1]
dataB6full <- eval(parse(text=paste0("dataB6full", t, condition)))
dataB6 <- eval(parse(text=paste0("dataB6", t, condition)))
dataBXDfull <- eval(parse(text=paste0("dataBXDfull", t, condition)))
dataBXD <- eval(parse(text=paste0("dataBXD", t, condition)))

# check compability of number of genes (rows)
ifelse(nrow(dataB6)==nrow(dataB6full) & nrow(dataBXD)==nrow(dataBXDfull), "Same number of genes", "Different number of genes!")

# join FastQTL statistics and q-values
eQTLB6 <- cbind(dataB6full, dataB6)
eQTLBXD <- cbind(dataBXDfull, dataBXD)

# define colnames
colnames(eQTLB6) <- c("nb_variants_tested", "MLE1_beta", "MLE2_beta", "TBD", "IDvariant", "distance", "pvalue", "slope", "pvalue_firstpermutation", "pvalue_secondpermutation", "adjustedpvalue", "marker")
colnames(eQTLBXD) <- c("nb_variants_tested", "MLE1_beta", "MLE2_beta", "TBD", "IDvariant", "distance", "pvalue", "slope", "pvalue_firstpermutation", "pvalue_secondpermutation", "adjustedpvalue", "marker")


## @knitr quantify

# find common genes, not na with any reference
commongenes <- intersect(rownames(eQTLB6)[!is.na(eQTLB6[, "slope"])], rownames(eQTLBXD)[!is.na(eQTLBXD[, "slope"])])
length(commongenes)

a <- eQTLB6[commongenes, "adjustedpvalue"]
b <- eQTLBXD[commongenes, "adjustedpvalue"]
diff <- abs((log10(a)-log10(b)))
plot(-log10(a), -log10(b), pch=21, bg=rgb(1,0,0, diff/max(diff)), main=paste(tissue, condition), xlab="-log10(qvalue) with B6 reference", ylab="-log10(qvalue) with BXD reference", las=1)
abline(a=0, b=1, lty=1)
abline(h=-log10(0.05), lty=2, col="darkred")
abline(v=-log10(0.05), lty=2, col="darkred")
plot(-log10(a), -log10(b), pch=21, bg=rgb(1,0,0, diff/max(diff)), main=paste(tissue, condition), xlab="-log10(qvalue) with B6 reference", ylab="-log10(qvalue) with BXD reference", las=1)
text(-log10(a), -log10(b), labels=commongenes, cex=0.7, pos=3)
abline(a=0, b=1, lty=1)
abline(h=-log10(0.05), lty=2, col="darkred")
abline(v=-log10(0.05), lty=2, col="darkred")

# plot slope comparison
a <- eQTLB6[commongenes, "slope"]
b <- eQTLBXD[commongenes, "slope"]
diff <- abs(a-b)
plot(a, b, pch=21, bg=rgb(1,0,0, diff/max(diff)), main=paste(tissue, condition), xlab="slope with B6 reference", ylab="slope with BXD reference", las=1)
abline(a=0, b=1, lty=1)
abline(h=0, lty=2)
abline(v=0, lty=2)
plot(a, b, pch=21, bg=rgb(1,0,0, diff/max(diff)), main=paste(tissue, condition), xlab="slope with B6 reference", ylab="slope with BXD reference", las=1)
text(a, b, labels=commongenes, cex=0.7, pos=3)
abline(a=0, b=1, lty=1)
abline(h=0, lty=2)
abline(v=0, lty=2)

# How many significant eQTLs are in common between the 2 references ?
sig_eQTLB6 <- commongenes[which(eQTLB6[commongenes, "adjustedpvalue"]<0.05)]
sig_eQTLBXD <- commongenes[which(eQTLBXD[commongenes, "adjustedpvalue"]<0.05)]
length(sig_eQTLB6)
length(sig_eQTLBXD)
mean(c(length(sig_eQTLB6), length(sig_eQTLBXD)))
commonsig <- intersect(sig_eQTLB6, sig_eQTLBXD)
length(commonsig)
# percentage of  significant eQTL in common?
length(commonsig) / mean(c(length(sig_eQTLB6), length(sig_eQTLBXD))) * 100

# how many genes have the same marker?
changeofmarker <- commongenes[eQTLB6[commongenes, "marker"]!=eQTLBXD[commongenes, "marker"]]
length(commongenes[eQTLB6[commongenes, "marker"]==eQTLBXD[commongenes, "marker"]])
length(commongenes)
# percentage of genes with the same marker?
length(commongenes[eQTLB6[commongenes, "marker"]==eQTLBXD[commongenes, "marker"]]) / length(commongenes) * 100

# how many genes are significant in both references?
length(setdiff(commonsig, changeofmarker))
length(setdiff(commonsig, changeofmarker)) / mean(c(length(sig_eQTLB6), length(sig_eQTLBXD))) * 100

##intersect(commonsig, changeofmarker)
##intersect(sig_eQTLB6, changeofmarker)
##intersect(sig_eQTLBXD, changeofmarker)

# Is this marginal differences (around FDR 0.05)
uniqsigB6 <- setdiff(sig_eQTLB6, sig_eQTLBXD)
length(uniqsigB6)
uniqsigBXD <- setdiff(sig_eQTLBXD, sig_eQTLB6)
length(uniqsigBXD)

hist(-log10(eQTLB6[uniqsigB6, "adjustedpvalue"]), breaks=50, las=1)
abline(v=-log10(0.05), lty=2, col="darkred")
hist(-log10(eQTLBXD[uniqsigBXD, "adjustedpvalue"]), breaks=50, las=1)
abline(v=-log10(0.05), lty=2, col="darkred")

# eQTL with same sign of slope
changeofsign <- commongenes[sign(eQTLB6[commongenes, "slope"])!=sign(eQTLBXD[commongenes, "slope"])]



##intersect(commonsig, changeofsign)
##intersect(commonsig, changeofmarker)
##intersect(changeofsign, changeofmarker)
##setdiff(setdiff(commonsig, changeofsign), changeofmarker)

# significant eQTLs maintained
##setdiff(setdiff(commonsig, changeofsign), changeofmarker)
length(setdiff(setdiff(commonsig, changeofsign), changeofmarker))
length(commonsig)
length(setdiff(setdiff(commonsig, changeofsign), changeofmarker)) / length(commonsig) * 100
# eQTLs maintained
##setdiff(setdiff(commongenes, changeofsign), changeofmarker)
length(setdiff(setdiff(commongenes, changeofsign), changeofmarker))
length(commongenes)
length(setdiff(setdiff(commongenes, changeofsign), changeofmarker)) / length(commongenes) * 100


a <- eQTLB6[commongenes, "adjustedpvalue"]
b <- eQTLB6[commongenes, "slope"]
plot(-log10(a), b, pch=21, bg=rgb(0,0,0, 0.1), main=paste(tissue, condition), xlab="-log10(qvalue) with B6 reference", ylab="slope with B6 reference", las=1)
abline(h=0, lty=1)
abline(v=-log10(0.05), lty=2, col="darkred")
plot(-log10(a), b, pch=21, bg=rgb(0,0,0, 0.1), main=paste(tissue, condition), xlab="-log10(qvalue) with B6 reference", ylab="slope with B6 reference", las=1)
abline(h=0, lty=1)
abline(v=-log10(0.05), lty=2, col="darkred")
text(-log10(a), b, labels=commongenes, cex=0.7, pos=3)

a <- eQTLBXD[commongenes, "adjustedpvalue"]
b <- eQTLBXD[commongenes, "slope"]
plot(-log10(a), b, pch=21, bg=rgb(0,0,0, 0.1), main=paste(tissue, condition), xlab="-log10(qvalue) with BXD reference", ylab="slope with BXD reference", las=1)
abline(h=0, lty=1)
abline(v=-log10(0.05), lty=2, col="darkred")
plot(-log10(a), b, pch=21, bg=rgb(0,0,0, 0.1), main=paste(tissue, condition), xlab="-log10(qvalue) with BXD reference", ylab="slope with BXD reference", las=1)
abline(h=0, lty=1)
abline(v=-log10(0.05), lty=2, col="darkred")
text(-log10(a), b, labels=commongenes, cex=0.7, pos=3)

# plot diagnostic for p-values
B6 <- eQTLB6[, "pvalue"]
BXD <- eQTLBXD[, "pvalue"]
hist(B6, breaks=50, col=rgb(1,0,0,0.5), xlab="p-values", 
     ylab="frequency", main=paste("Diagnostic p-values plot", tissue, condition), las=1)
hist(BXD, breaks=50, col=rgb(0,0,1,0.5), add=TRUE)
# add legend
legend("topright", legend=c("B6","BXD"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15)

# plot diagnostic for q-values
B6 <- eQTLB6[, "adjustedpvalue"]
BXD <- eQTLBXD[, "adjustedpvalue"]
hist(B6, breaks=50, col=rgb(1,0,0,0.5), xlab="adjusted p-values", 
     ylab="frequency", main=paste("Diagnostic adjusted p-values plot", tissue, condition), las=1)
hist(BXD, breaks=50, col=rgb(0,0,1,0.5), add=TRUE)
# add legend
legend("topright", legend=c("B6","BXD"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15)

# percentage of significant eQTLs over expressed genes
length(which(eQTLB6[, "adjustedpvalue"]<0.05)) / length(eQTLB6[, "adjustedpvalue"]) * 100
length(which(eQTLBXD[, "adjustedpvalue"]<0.05)) / length(eQTLBXD[, "adjustedpvalue"]) * 100


## @knitr qvalue_package_diagnostic

qB6 <- qvalue::qvalue(eQTLB6$pvalue_secondpermutation)
qBXD <- qvalue::qvalue(eQTLBXD$pvalue_secondpermutation)

summary(qB6)
hist(qB6)
plot(qB6)

plot(qB6$pvalues, qB6$qvalues)
abline(a=0, b=1, lty=2)

summary(qBXD)
hist(qBXD)
plot(qBXD)

plot(qBXD$pvalues, qBXD$qvalues)
abline(a=0, b=1, lty=2)

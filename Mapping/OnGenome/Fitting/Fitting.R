setwd("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/Fitting")

# try to fit y (=uniquely mapped %) according to x (= number mismatches allowed)
# GOAL to do some categories, spot the extremes


# load stats for all samples
B6mm10all <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenome/MappingStatisticsB6mm10_primaryassembly.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")
D2all <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenome/MappingStatisticsD2.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")


# test chromosomes and scaffolds to include
##B6mm10primaryassembly <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenome/stats_test/MappingStatisticsB6mm10primaryassembly_B6D2.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")
##B6mm10primaryassemblynoscaffolds <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenome/stats_test/MappingStatisticsB6mm10primaryassemblynoscaffolds_B6D2.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")
##B6mm10majorchromosomes <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenome/stats_test/MappingStatisticsB6mm10majorchromosomes_B6D2.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")
##D2majorchromosomes <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenome/stats_test/MappingStatisticsD2majorchromosomes_B6D2.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")

# plotting all in 1 graph
plot("","", las=1, bty="l", pch=19, xlim=range(0:10), ylim=c(floor(min(B6mm10all,D2all)),ceiling(max(B6mm10all,D2all))), xlab="# mismatch allowed", ylab="Uniquely mapped %")
##par(mfrow=c(2,5))
for(mysample in rownames(B6mm10all)){
  ##print(mysample)
  # on B6 mm10 reference
  points(c(0:5,10),B6mm10all[mysample,], las=1, bty="l", pch=19)
  # on D2 reference
  points(c(0:5,10),D2all[mysample,], las=1, bty="l", pch=19, col="tan")
}


B6mm10all_m <- t(B6mm10all)
D2all_m <- t(D2all)

# quick linear model
x <- as.numeric(gsub("Uniq", "", rownames(B6mm10all_m)))
lm_1 <- lm(B6mm10all_m[,1] ~ x)
summary(lm_1)
# plot fitting
plot(B6mm10all_m[,1] ~ x, ylim=c(0,80), pch=19)
lines(x, predict(lm_1, data.frame(x)), col='darkblue', lwd=1)
# check residuals
plot(lm_1)

# second degree
lm_2 <- lm(B6mm10all_m[,168] ~ x+I(x^2))
summary(lm_2)
# plot fitting
plot(B6mm10all_m[,168] ~ x, ylim=c(0,80), pch=19)
lines(x, predict(lm_2, data.frame(x)), col='darkblue', lwd=1)
# check residuals
plot(lm_2)

lm_2D <- lm(D2all_m[,168] ~ x+I(x^2))
summary(lm_2)
# plot fitting
plot(D2all_m[,168] ~ x, ylim=c(0,80), pch=19)
lines(x, predict(lm_2D, data.frame(x)), col='darkblue', lwd=1)
# check residuals
plot(lm_2D)

# 3rd degree
lm_3 <- lm(B6mm10all_m[,1] ~ x+I(x^2)+I(x^3))
summary(lm_3)
# plot fitting
plot(B6mm10all_m[,1] ~ x, ylim=c(0,80), pch=19)
lines(x, predict(lm_3, data.frame(x)), col='darkblue', lwd=1)
# check residuals
plot(lm_3)

# Conclusion: second degree is kind of ok.

# create a list to register lm for the different samples

# all at once
x <- as.numeric(gsub("Uniq", "", rownames(B6mm10all_m)))
lm_B6mm10_all <- lm(B6mm10all_m ~ x+I(x^2))
##lm_B6mm10_all <- lm(B6mm10all_m ~ x+I(x^2)+I(x^3))
lmsB6 <- summary(lm_B6mm10_all)
lm_D2_all <- lm(D2all_m ~ x+I(x^2))
##lm_D2_all <- lm(D2all_m ~ x+I(x^2)+I(x^3))
lmsD2 <- summary(lm_D2_all)

par(mfrow=c(2,5))
for(mysample in rownames(B6mm10all)){
  plot("","", las=1, bty="l", pch=19, xlim=range(x), ylim=c(floor(min(B6mm10all,D2all)),ceiling(max(B6mm10all,D2all))), xlab="# mismatch allowed", ylab="Uniquely mapped %", main=mysample)
  # on B6 mm10 reference
  points(x,B6mm10all[mysample,], las=1, bty="l", pch=19)
  # on D2 reference
  points(x,D2all[mysample,], las=1, bty="l", pch=19, col="tan")
  # plot B6 fitting
  lines(x, lm_B6mm10_all$fitted.values[,mysample], col='darkblue', lwd=1)
  # plot D2 fitting
  lines(x, lm_D2_all$fitted.values[,mysample], col='chocolate', lwd=1)
}

B6minusD2all <- B6mm10all_m-D2all_m
##lm_B6minusD2_all <- lm(B6minusD2all ~ x+I(x^2)+I(x^3))
lm_B6minusD2_all <- lm(B6minusD2all ~ x+I(x^2))
lmsB6minusD2 <- summary(lm_B6minusD2_all)

##par(mfrow=c(1,1))
##plot("","", las=1, ylim=c(-10,10), xlim=range(x), bty="l", pch=19)
for(mysample in rownames(B6mm10all)){
  plot("","", las=1, ylim=c(-10,10), xlim=range(x), bty="l", pch=19, main=mysample)
  # plot reference
  points(x,B6minusD2all[,mysample], las=1, bty="l", pch=19)
  # plot fitting
  lines(x, lm_B6minusD2_all$fitted.values[,mysample], col='darkblue', lwd=1)
  # find 0
  print(c(mysample, solve(polynomial(lm_B6minusD2_all$coefficients[,mysample]))))
}


library(polynom)
solve(polynomial(lm_B6minusD2_all$coefficients[,mysample]))
# source: https://r.789695.n4.nabble.com/solving-cubic-quartic-equations-non-iteratively-td999092.html

####################################################################

par(mfrow=c(1,1))
##plot(B6minusD2all[1,], pch=19)
plot(B6minusD2all[4,], pch="")
text(B6minusD2all[4,], labels(B6minusD2all[4,]))
for (i in rownames(B6minusD2all)){
  plot(B6minusD2all[i,], pch="", main=i, ylim=c(-10,10))
  abline(h=0, col="darkred")
  text(B6minusD2all[i,], labels(B6minusD2all[i,]))
}

parental <- grep("B[612]", colnames(B6minusD2all))
for (i in rownames(B6minusD2all)){
  plot(B6minusD2all[i,parental], pch="", main=i, ylim=c(-10,10))
  abline(h=0, col="darkred")
  text(B6minusD2all[i,parental], labels(B6minusD2all[i,parental]))
  print(c(i, mean(B6minusD2all[i,parental])))
}

# Uniq4 seems to present the least bias for a reference
plot(B6minusD2all[4,], pch="")
text(B6minusD2all[4,], labels(B6minusD2all[4,]))


#######################################

B6minusD2_MC <- t(B6mm10majorchromosomes-D2majorchromosomes)
for (i in rownames(B6minusD2_MC)){
  plot(B6minusD2_MC[i,], pch="", main=i, ylim=c(-10,10))
  abline(h=0, col="darkred")
  text(B6minusD2_MC[i,], labels(B6minusD2_MC[i,]))
  print(c(i, mean(B6minusD2_MC[i,])))
}

B6minusD2_pa <- t(B6mm10primaryassembly-D2all[rownames(B6mm10primaryassembly),])
for (i in rownames(B6minusD2_pa)){
  plot(B6minusD2_pa[i,], pch="", main=i, ylim=c(-10,10))
  abline(h=0, col="darkred")
  text(B6minusD2_pa[i,], labels(B6minusD2_pa[i,]))
  print(c(i, mean(B6minusD2_pa[i,])))
}

B6minusD2_pans <- t(B6mm10primaryassemblynoscaffolds-D2all[rownames(B6mm10primaryassemblynoscaffolds),])
for (i in rownames(B6minusD2_pans)){
  plot(B6minusD2_pans[i,], pch="", main=i, ylim=c(-10,10))
  abline(h=0, col="darkred")
  text(B6minusD2_pans[i,], labels(B6minusD2_pans[i,]))
  print(c(i, mean(B6minusD2_pans[i,])))
}

####################################################################

plot(lm_B6mm10_all$coefficients[1,], lm_B6mm10_all$coefficients[2,], pch="")
text(lm_B6mm10_all$coefficients[1,], lm_B6mm10_all$coefficients[2,], labels(lm_B6mm10_all$coefficients[1,]))
plot(lm_B6mm10_all$coefficients[1,], lm_B6mm10_all$coefficients[3,], pch="")
text(lm_B6mm10_all$coefficients[1,], lm_B6mm10_all$coefficients[3,], labels(lm_B6mm10_all$coefficients[1,]))
plot(lm_B6mm10_all$coefficients[2,], lm_B6mm10_all$coefficients[3,], pch="")
text(lm_B6mm10_all$coefficients[2,], lm_B6mm10_all$coefficients[3,], labels(lm_B6mm10_all$coefficients[2,]))

# check extreme values
##lm_all$coefficients[,c("32nsd","32sd","B61nsd","B61sd","B62nsd","B62sd", "BDnsd","BDsd")]
##lm_all$coefficients[,c("L29tnsd","L29tsd","LB61nsd","LB61sd","LB62nsd","LB62sd", "LDB2sd","LDB1nsd")]

for(i in 1:3){
  print(labels(c(which.min(lm_all$coefficients[i,]),which.max(lm_all$coefficients[i,]))))
}
lm_all$coefficients[,c("LDB2sd", "B62nsd", "LB62sd", "32sd", "32sd", "L29tnsd")]

# in cortex
for(i in 1:3){
  print(labels(c(which.min(lm_all$coefficients[i,1:88]),which.max(lm_all$coefficients[i,1:88]))))
}
lm_all$coefficients[,c("32nsd", "B62nsd", "B61nsd", "32sd", "32sd", "B61nsd")]

# in liver
for(i in 1:3){
  print(labels(c(which.min(lm_all$coefficients[i,89:176]),which.max(lm_all$coefficients[i,89:176]))))
}
lm_all$coefficients[,c("LDB2sd", "LB62sd","LB62sd", "LDB1nsd", "LDB1nsd", "L29tnsd")]

plot("","", las=1, bty="l", pch=19, xlim=range(0:5), ylim=c(floor(min(B6mm10all,D2all)),ceiling(max(B6mm10all,D2all))), xlab="# mismatch allowed", ylab="Uniquely mapped %", main=mysample)
##par(mfrow=c(2,5))
for(mysample in rownames(B6mm10all)){
  # on B6 mm10 reference
  points(c(0:5),B6mm10all[mysample,], las=1, bty="l", pch=19)
  # on D2 reference
  points(c(0:5),D2all[mysample,], las=1, bty="l", pch=19, col="tan")
}
plot(rep(0:5, length(c("LDB2sd", "LB62sd", "LB62sd", "LDB1nsd", "LDB1nsd", "L29tnsd"))), B6mm10all_m[, c("LDB2sd", "LB62sd","LB62sd", "LDB1nsd", "LDB1nsd", "L29tnsd")], pch=19, ylim=c(0,70))
points(rep(0:5, length(c("LDB2sd", "LB62sd", "LB62sd", "LDB1nsd", "LDB1nsd", "L29tnsd"))), predict(lm_all, data.frame(x)), col='blue', lwd=2)

# get adjusted R^2
R2 <- vector(length=ncol(B6mm10all_m), mode="numeric")
for(mysampleidx in 1:ncol(B6mm10all_m)){
  R2[mysampleidx] <- as.numeric(lms[[mysampleidx]]["adj.r.squared"])
}
R2



##############################################################

# second degree
x <- as.numeric(gsub("Uniq", "", rownames(D2all_m)))
lm_D2_all <- lm(D2all_m ~ x+I(x^2))
summary(lm_D2_all)

# in cortex
for(i in 1:3){
  print(labels(c(which.min(lm_D2_all$coefficients[i,1:88]),which.max(lm_D2_all$coefficients[i,1:88]))))
}
lm_D2_all$coefficients[,c("B61nsd", "DB2sd", "DB1nsd", "43sd", "43sd", "DB1nsd")]

# in liver
for(i in 1:3){
  print(labels(c(which.min(lm_D2_all$coefficients[i,89:176]),which.max(lm_D2_all$coefficients[i,89:176]))))
}
lm_D2_all$coefficients[,c("LB61nsd", "L81sd","LDB2sd", "L73nsd", "L73nsd", "L44sd")]

###############################################################

# second degree
x <- as.numeric(gsub("Uniq", "", rownames(D2all_m)))
y <- B6mm10all_m - D2all_m
lm_diff_all <- lm(y ~ x+I(x^2))
summary(lm_diff_all)

# in cortex
for(i in 1:3){
  print(labels(c(which.min(lm_diff_all$coefficients[i,1:88]),which.max(lm_diff_all$coefficients[i,1:88]))))
}
lm_diff_all$coefficients[,c("DB2sd", "B61sd", "B61sd", "DB2sd", "DB1nsd", "B62nsd")]

# in liver
for(i in 1:3){
  print(labels(c(which.min(lm_diff_all$coefficients[i,89:176]),which.max(lm_diff_all$coefficients[i,89:176]))))
}
lm_diff_all$coefficients[,c("LDB2sd", "LB62nsd", "LB62nsd", "LDB2nsd", "LDB2nsd", "LB62nsd")]

plot(lm_diff_all$coefficients[1,], lm_diff_all$coefficients[2,], pch="")
text(lm_diff_all$coefficients[1,], lm_diff_all$coefficients[2,], labels(lm_diff_all$coefficients[1,]))
plot(lm_diff_all$coefficients[1,], lm_diff_all$coefficients[3,], pch="")
text(lm_diff_all$coefficients[1,], lm_diff_all$coefficients[3,], labels(lm_diff_all$coefficients[1,]))
plot(lm_diff_all$coefficients[2,], lm_diff_all$coefficients[3,], pch="")
text(lm_diff_all$coefficients[2,], lm_diff_all$coefficients[3,], labels(lm_diff_all$coefficients[2,]))

# plotting all in 1 graph
plot("","", las=1, bty="l", pch=19, xlim=range(x), ylim=c(floor(min(y)),ceiling(max(y))), xlab="# mismatch allowed", ylab="Uniquely mapped %")
##par(mfrow=c(2,5))
for(mysample in colnames(y)){
  # on B6 mm10 - D2
  points(x,y[,mysample], las=1, bty="l", pch=19)
}

# plotting all in 1 graph
plot("","", las=1, bty="l", pch=19, xlim=range(x), ylim=c(floor(min(y)),ceiling(max(y))), xlab="# mismatch allowed", ylab="Uniquely mapped %")
##par(mfrow=c(2,5))
for(mysample in c("LDB2sd", "LB62nsd", "LB62nsd", "LDB2nsd", "LDB2nsd", "LB62nsd")){
  # on B6 mm10 - D2
  points(x,y[,mysample], las=1, bty="l", pch=19)
}



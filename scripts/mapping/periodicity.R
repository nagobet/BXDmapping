# Checking periodicity
# Plotting + testing liver uniquely mapped values for STAR default on genome + transcriptome annotation
# 12 February 2020

# How to fit a sine curve?
# example from source: https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
y <- c(11.622967, 12.006081, 11.760928, 12.246830, 12.052126, 12.346154, 12.039262, 12.362163, 12.009269, 11.260743, 10.950483, 10.522091,  9.346292,  7.014578,  6.981853,  7.197708,  7.035624,  6.785289, 7.134426,  8.338514,  8.723832, 10.276473, 10.602792, 11.031908, 11.364901, 11.687638, 11.947783, 12.228909, 11.918379, 12.343574, 12.046851, 12.316508, 12.147746, 12.136446, 11.744371,  8.317413, 8.790837, 10.139807,  7.019035,  7.541484,  7.199672,  9.090377,  7.532161,  8.156842,  9.329572, 9.991522, 10.036448, 10.797905)
t <- 18:65
ssp <- spectrum(y)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
reslm <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t))
summary(reslm)

rg <- diff(range(y))
plot(y~t,ylim=c(min(y)-0.1*rg,max(y)+0.1*rg))
lines(fitted(reslm)~t,col=4,lty=2)   # dashed blue line is sin fit

# including 2nd harmonic really improves the fit
reslm2 <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t)+sin(4*pi/per*t)+cos(4*pi/per*t))
summary(reslm2)
lines(fitted(reslm2)~t,col=3)    # solid green line is periodic with second harmonic
################################

# load my data
map_B6mm10 <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/2-pass/MappingStatisticsB6mm10.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")
map_B6mm9 <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/2-pass/MappingStatisticsB6mm9.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")
map_D2 <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/2-pass/MappingStatisticsD2.tsv", header=TRUE, stringsAsFactors=FALSE, row.names="Sample_Name")

sratable <- read.table("F:/BXD/data/transcriptome/RNAseq/SraRunTable.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
uniqSratmp <- unique(sratable[, c("Sample_Name", "condition", "genotype", "tissue")])
# load table to convert mouse lines names
MouseLines <- read.table("F:/BXD/data/ConvertLineNames.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
# clean mouse lines names in sra table
uniqSra <- uniqSratmp
rownames(uniqSra) <- 1:nrow(uniqSra)
uniqSra$genotype <- MouseLines$V1[match(uniqSratmp$genotype,MouseLines$V2)]
rm(uniqSratmp, MouseLines)
# Apply it to my data

##y <- B6mm10all$Uniq0[(176/2+1):176]
y <- map_B6mm10$Uniq[(172/2+1):172]
##t <- (176/2+1):176
t <- (172/2+1):172

ssp <- spectrum(y)
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
reslm <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t))
##reslm <- lm(y ~ sin(2*pi/per*t))
summary(reslm)

rg <- diff(range(y))
plot(y~t,ylim=c(min(y)-0.1*rg,max(y)+0.1*rg), pch=19, xaxt='n')
lines(fitted(reslm)~t,lty=2)   # dashed blue line is sin fit
axis(side=1, at=(172/2+1):172, labels=uniqSra$genotype[match(rownames(map_B6mm10)[(172/2+1):172],uniqSra$Sample_Name)], las=2)


# 24 February 2020
# plot separately conditions
par(mfrow=c(2,2))
for(c in c("Control", "Sleep Deprived")){
  for(t in c("Liver", "Cortex")){
    myselection <- uniqSra$condition==c&uniqSra$tissue==t
    myorder <- order(uniqSra$genotype[match(rownames(map_B6mm10)[myselection],uniqSra$Sample_Name)])
    mylabels <- uniqSra$genotype[match(rownames(map_B6mm10)[myselection],uniqSra$Sample_Name)]
    mylabels <- mylabels[myorder]
    myvalues <- map_B6mm10$Uniq[myselection]
    myvalues <- myvalues[myorder]
    plot(1:43, myvalues, pch=19, xlab="", xaxt='n', las=2, main=paste(c, t))
    axis(side=1, at=1:43, labels=mylabels, las=2, cex.axis=0.7)
    ssp <- spectrum(myvalues, plot=FALSE)
    per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
    reslm <- lm(myvalues ~ sin(2*pi/per*1:43)+cos(2*pi/per*1:43))
    sreslm <- summary(reslm)
    pval <- pf(sreslm$fstatistic[1],sreslm$fstatistic[2],sreslm$fstatistic[3], lower.tail=FALSE)
    mtext(paste("p-value:", format(pval, digits=4, scientific=TRUE)),side=3, adj=1)
    lines(1:43, fitted(reslm), lty=2)   # dashed line is sin fit
    print(summary(reslm))
  }
}

# random lines order
for(c in c("Control", "Sleep Deprived")){
  for(t in c("Liver", "Cortex")){
    myselection <- uniqSra$condition==c&uniqSra$tissue==t
    myorder <- sample(1:length(uniqSra$genotype[match(rownames(map_B6mm10)[myselection],uniqSra$Sample_Name)]))
    mylabels <- uniqSra$genotype[match(rownames(map_B6mm10)[myselection],uniqSra$Sample_Name)]
    mylabels <- mylabels[myorder]
    myvalues <- map_B6mm10$Uniq[myselection]
    myvalues <- myvalues[myorder]
    plot(1:43, myvalues, pch=19, xlab="", xaxt='n', las=2, main=paste(c, t))
    axis(side=1, at=1:43, labels=mylabels, las=2, cex.axis=0.7)
    ssp <- spectrum(myvalues, plot=FALSE)
    per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
    reslm <- lm(myvalues ~ sin(2*pi/per*1:43)+cos(2*pi/per*1:43))
    sreslm <- summary(reslm)
    pval <- pf(sreslm$fstatistic[1],sreslm$fstatistic[2],sreslm$fstatistic[3], lower.tail=FALSE)
    mtext(paste("p-value:", format(pval, digits=4, scientific=TRUE)),side=3, adj=1)
    lines(1:43, fitted(reslm), lty=2)   # dashed line is sin fit
    print(summary(reslm))
  }
}
# order mouse lines by uniquely mapped percentage
for(c in c("Control", "Sleep Deprived")){
  for(t in c("Liver", "Cortex")){
    # filter for tissue and condition
    myselection <- uniqSra$condition==c&uniqSra$tissue==t
    mylabels <- uniqSra$genotype[match(rownames(map_B6mm10)[myselection],uniqSra$Sample_Name)]
    myvalues <- map_B6mm10$Uniq[myselection]
    # order mouse lines
    myorder <- order(myvalues)
    mylabels <- mylabels[myorder]
    myvalues <- myvalues[myorder]
    # plot
    plot(1:43, myvalues, pch=19, xlab="", xaxt='n', las=2, main=paste(c, t))
    axis(side=1, at=1:43, labels=mylabels, las=2, cex.axis=0.7)
  }
}

# plot Uniq values by sample all conditions and tissues on one plot
par(mfrow=c(1,1))
plot(c(1,43), c(60,95), pch="", las=2, xaxt="n",xlab="", ylab="Uniquely mapped (%)")
axis(side=1, at=1:43, labels=uniqSra$genotype[match(rownames(map_B6mm10)[myselection],uniqSra$Sample_Name)], las=2, cex.axis=0.9)
points(rep(sort(rep(1:43,2)),2),map_B6mm10$Uniq, pch=20, col=c(rep(1:2,43),rep(3:4,43)))
legend("bottomleft",legend=c("NSD Cortex", "SD Cortex","NSD Liver", "SD Liver"), pch=20, col=1:4)

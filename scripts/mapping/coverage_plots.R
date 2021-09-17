#plot coverage
##read.table("F:/BXD/LDB1exacttmp.cov", header=FALSE)
##cov <- read.table("F:/BXD/tmp5000", header=FALSE)
##cov <- read.table("F:/BXD/tmp100000", header=FALSE)

##covdataB6 <- read.table("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/B6mm10_majorchromosomes/LB61nsd_exactunique.cov", header=FALSE)
##covdataD2 <- read.table("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/D2_majorchromosomes/LB61nsd_exactunique.cov", header=FALSE)

##plot(cov$V2,cov$V3, type="line")
##plot(cov$V2)

##plot(cov_2$V2,log10(cov_2$V3), type="line")

##png(filename=paste("chromosome_", chr, ".png", sep=""), width=480, height=190)
##plot(1:10, xlab="x-label", main="My nice title")
##dev.off()

# function to produce plot
plotCov <- function(sample="LB61nsd", refB6="B6mm10_primaryassembly", refD2="D2"){
##chr <- 1
##chr <- 2
  for(chr in c(1:19,"X")){
  ##for(chr in c(1,6,8)){
    print(paste("loading data for chr",chr))
    covdataB6 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refB6, "/", sample, "_", chr,".cov", sep=""), header=FALSE)
    covdataD2 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refD2, "/", sample, "_", chr,".cov", sep=""), header=FALSE)
    ##covdataB6 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refB6, "/", sample, "_", chr,".cov", sep=""), header=FALSE)
    ##covdataD2 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refD2, "/", sample, "_", chr,".cov", sep=""), header=FALSE)
    ##png(filename=paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", sample, "_", chr, ".png", sep=""), width=900, height=300)
    print("data loaded, plotting")
    ##shift <- covdataB6$V2[covdataB6$V3==max(covdataB6$V3)]-covdataD2$V2[covdataD2$V3==max(covdataD2$V3)]
    ##shift <- 0
    shift <- max(covdataB6$V2)-max(covdataD2$V2)
    png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/Coverage/plots_logE_shiftEnd/", sample, "_", chr, ".png", sep=""), width=900, height=300)
    plot(x="", y="", xlim=c(0,max(covdataB6$V2,covdataD2$V2)), ylim=range(-11,11), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log(coverage)")
    lines(covdataB6$V2,log(covdataB6$V3), col="black")
    ##lines(covdataD2$V2,-log(covdataD2$V3), col="tan3")
    lines(covdataD2$V2+shift,-log(covdataD2$V3), col="tan3")
    text(1,c(9,-9), labels=c("B6","D2"), col=c("black","tan3"))
    dev.off()
  }
}

##plot(covdataB6$V2[covdataB6$V1==chr],covdataB6$V3[covdataB6$V1==chr], ylim=range(-covdataD2$V3[covdataD2$V1==chr],covdataB6$V3[covdataB6$V1==chr]), type="line", main=paste("chromosome", chr), xlab="position [bp]", ylab="log10(coverage)")
##lines(covdataD2$V2[covdataD2$V1==chr],-covdataD2$V3[covdataD2$V1==chr], col="tan4")

plotCov(sample="LB61nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="LB62nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="LDB1nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="LDB2nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="B61nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="B62nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="DB1nsd", refB6="B6mm10_primaryassembly", refD2="D2")
plotCov(sample="DB2nsd", refB6="B6mm10_primaryassembly", refD2="D2")

##png(filename=paste("/mnt/nas/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/test", sep=""), width=900, height=300)
##plot(1:10)
##dev.off()

# trying to catch the shift
sample="LDB1nsd"
##sample <- "B61nsd"
refB6="B6mm10_primaryassembly"
refD2="D2"
chr <- 14
##chr <- 1
covdataB6 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refB6, "/", sample, "_", chr,".cov", sep=""), header=FALSE)
covdataD2 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refD2, "/", sample, "_", chr,".cov", sep=""), header=FALSE)

##summary(covdataB6$V3)
##summary(covdataD2$V3)

##covdataB6$V2[covdataB6$V3==4191]

shift <- covdataB6$V2[covdataB6$V3==max(covdataB6$V3)]-covdataD2$V2[covdataD2$V3==max(covdataD2$V3)]

##min(covdataB6$V2)-min(covdataD2$V2)
##shift <- max(covdataB6$V2)-max(covdataD2$V2)

##shift <- mean(c(min(covdataB6$V2)-min(covdataD2$V2),max(covdataB6$V2)-max(covdataD2$V2)))


plot(x="", y="", xlim=c(0,max(covdataB6$V2,covdataD2$V2)), ylim=range(-11,11), main=paste(sample, " chromosome", chr), xlab="position [bp]", ylab="log10(coverage)")
lines(covdataB6$V2,log(covdataB6$V3), col="black")
lines(covdataD2$V2+shift,-log(covdataD2$V3), col="tan3")

# try with other transformation than log
# sqrt
plot(x="", y="", xlim=c(0,max(covdataB6$V2,covdataD2$V2)), ylim=range(-65,65), main=paste(sample, " chromosome", chr), xlab="position [bp]", ylab="coverage^0.5")
lines(covdataB6$V2,sqrt(covdataB6$V3), col="black")
lines(covdataD2$V2+shift,-sqrt(covdataD2$V3), col="tan3")

# no transformation
plot(x="", y="", xlim=c(0,max(covdataB6$V2,covdataD2$V2)), ylim=range(-4500,4500), main=paste(sample, " chromosome", chr), xlab="position [bp]", ylab="coverage")
lines(covdataB6$V2,covdataB6$V3, col="black")
lines(covdataD2$V2+shift,-covdataD2$V3, col="tan3")

################################

# plot coverage and annotations

chr <- 1
sample="LB61nsd"
refD2="D2"
covdataD2 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", refD2, "/", sample, "_", chr,".cov", sep=""), header=FALSE)
covdataB6mm10 <- read.table(paste("F:/BXD/data/transcriptome/3_coverage_bedtools/OnGenome/", "B6mm10", "/", sample, "_", chr,".cov", sep=""), header=FALSE)

shift <- 0
geneD2 <- read.table("F:/BXD/references/geneD2.tab", header=FALSE, sep="\t", stringsAsFactors=FALSE)
geneB6mm10 <- read.table("F:/BXD/references/geneB6mm10.tab", header=FALSE, sep="\t", stringsAsFactors=FALSE)
exonD2 <- read.table("F:/BXD/references/exonD2.tab", header=FALSE, sep="\t", stringsAsFactors=FALSE)
exonB6mm10 <- read.table("F:/BXD/references/exonB6mm10.tab", header=FALSE, sep="\t", stringsAsFactors=FALSE)



png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/CovAnn_", sample, "_", chr, ".png", sep=""), width=900, height=300)
##plot(covdataD2$V2[1102335:2204670],log(covdataD2$V3[1102335:2204670]), col="tan3", ylim=range(-10,10), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log(coverage)", type="l")
plot(covdataD2$V2,log(covdataD2$V3), col="tan3", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataD2$V2))), ylim=range(-10,10), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log(coverage)", type="l")
rect(xleft=exonD2$V4[exonD2$V1==chr], ybottom=-5, xright=exonD2$V5[exonD2$V1==chr], ytop=-3, col="tan", border="tan")
rect(xleft=exonB6mm10$V4[exonB6mm10$V1==chr], ybottom=-9, xright=exonB6mm10$V5[exonB6mm10$V1==chr], ytop=-7, col="black", border="black")
rect(xleft=geneD2$V4[geneD2$V1==chr], ybottom=-4, xright=geneD2$V5[geneD2$V1==chr], ytop=-4, col="tan", border="tan")
rect(xleft=geneB6mm10$V4[geneB6mm10$V1==chr], ybottom=-8, xright=geneB6mm10$V5[geneB6mm10$V1==chr], ytop=-8, col="black", border="black")
dev.off()

# test transformation
# no transformation
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_D2_notransform.png", sep=""), width=900, height=300)
plot(covdataD2$V2,covdataD2$V3, col="tan3", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataD2$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="coverage", type="l")
dev.off()
# natural logarithm
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_D2_lognatural.png", sep=""), width=900, height=300)
plot(covdataD2$V2,log(covdataD2$V3), col="tan3", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataD2$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log_e(coverage)", type="l")
dev.off()
# log 10
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_D2_log10.png", sep=""), width=900, height=300)
plot(covdataD2$V2,log10(covdataD2$V3), col="tan3", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataD2$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log10(coverage)", type="l")
dev.off()
# log 2
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_D2_log2.png", sep=""), width=900, height=300)
plot(covdataD2$V2,log2(covdataD2$V3), col="tan3", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataD2$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log2(coverage)", type="l")
dev.off()
# square root
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_D2_sqrt.png", sep=""), width=900, height=300)
plot(covdataD2$V2,sqrt(covdataD2$V3), col="tan3", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataD2$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="sqrt(coverage)", type="l")
dev.off()

# test transformation
# no transformation
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_notransform.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2,covdataB6mm10$V3, col="black", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="coverage", type="l")
dev.off()
# natural logarithm
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_lognatural.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2,log(covdataB6mm10$V3), col="black", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log_e(coverage)", type="l")
dev.off()
# log 10
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_log10.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2,log10(covdataB6mm10$V3), col="black", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log10(coverage)", type="l")
dev.off()
# log 2
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_log2.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2,log2(covdataB6mm10$V3), col="black", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log2(coverage)", type="l")
dev.off()
# square root
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_sqrt.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2,sqrt(covdataB6mm10$V3), col="black", xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="sqrt(coverage)", type="l")
dev.off()

# filtering for low coverage
# justificable threshold
hist(covdataB6mm10$V3)
hist(log(covdataB6mm10$V3))
hist(sqrt(covdataB6mm10$V3))

# natural logarithm
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_lognatural_filter10.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2[log(covdataB6mm10$V3)>log(10)],log(covdataB6mm10$V3)[log(covdataB6mm10$V3)>log(10)], col="black", ylim=c(0,11), xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log_e(coverage)", type="l")
dev.off()
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_lognatural_filter100.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2[log(covdataB6mm10$V3)>log(100)],log(covdataB6mm10$V3)[log(covdataB6mm10$V3)>log(100)], col="black", ylim=c(0,11), xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="log_e(coverage)", type="l")
dev.off()

# square root
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_sqrt_filter10.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2[sqrt(covdataB6mm10$V3)>sqrt(10)],sqrt(covdataB6mm10$V3)[sqrt(covdataB6mm10$V3)>sqrt(10)], col="black", ylim=c(0,220), xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="sqrt(coverage)", type="l")
dev.off()
png(filename=paste("F:/BXD/analyses/Comparison_mapping_B6vsD2/OnGenome/coverage/", sample, "_", chr, "_on_B6mm10_sqrt_filter100.png", sep=""), width=900, height=300)
plot(covdataB6mm10$V2[sqrt(covdataB6mm10$V3)>sqrt(100)],sqrt(covdataB6mm10$V3)[sqrt(covdataB6mm10$V3)>sqrt(100)], col="black", ylim=c(0,220), xlim=c(0, max(c(geneB6mm10$V5, geneD2$V5, covdataB6mm10$V2))), main=paste("sample:", sample, " chromosome", chr), xlab="position [bp]", ylab="sqrt(coverage)", type="l")
dev.off()

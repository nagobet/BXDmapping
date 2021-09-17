# 15 November 2019

# GOAL: Compare kallisto mapping percentage in liver vs in cortex

RNA <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/GSE114845-expression.txt.txt", header=TRUE, stringsAsFactors=FALSE)
Meta <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/GSE114845-metadata.txt.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(Meta) <- Meta$Sample_geo_accession
rownames(RNA) <- RNA$gene_symbol

##countsCortexKallisto <- RNA[,Meta$Sample_geo_accession[Meta$tissue == "Cortex" ]]
##Meta_cortex <- Meta[Meta$tissue == "Cortex" ,]

##countsLiverKallisto <- RNA[,Meta$Sample_geo_accession[Meta$tissue == "Liver"]]
##Meta_liver <- Meta[Meta$tissue == "Liver",]


##head(countsCortexKallisto)
##apply(countsCortexKallisto, 2, sum)
##apply(countsLiverKallisto, 2, sum)


# load total sequences

##totseq <- read.table("F:/BXD/data/transcriptome/1_QC_fastqc/TotalSequences_allsamples.tab", header=TRUE, stringsAsFactors=FALSE)



# plot

# mapped percentage
##plot(c(apply(countsCortexKallisto, 2, sum),apply(countsLiverKallisto, 2, sum))/t(totseq)*100)

# mapped percentage cortex vs liver
##plot(apply(countsCortexKallisto, 2, sum),apply(countsLiverKallisto, 2, sum))
##abline(a=0,b=1, col="darkred")

# total sequences
##plot(t(totseq))

# total sequences cortex vs liver
##plot(t(totseq[,1:86]),t(totseq[,87:172]))
##abline(a=0,b=1, col="darkred")

##plot(t(totseq[,1:86])-t(totseq[,87:172]))


##apply(countsLiverKallisto, 2, sum)
##totseq[, 86+(60:70)]


#####################

# 4 December 2019

# load stats
n_pseudoaligned <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/n_pseudoaligned.tab", header=FALSE, stringsAsFactors=FALSE)
n_processed <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/n_processed.tab", header=FALSE, stringsAsFactors=FALSE)
n_unique <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/n_unique.tab", header=FALSE, stringsAsFactors=FALSE)

# quick visual comparison cortex vs liver
plot(n_unique$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Pseudomapping with Kallisto", xlab="Samples", ylab="Uniquely mapped percentage",
     las=1)
legend("bottomright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")
plot(n_pseudoaligned$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Pseudomapping with Kallisto", xlab="Samples", ylab="Pseudoaligned percentage",
     las=1)
legend("bottomright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")


# 13 December 2019
# plot the graphs with own-build indexes

# stats with D2 as a reference

# load stats
n_pseudoaligned <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/D2/n_pseudoaligned.tab", header=FALSE, stringsAsFactors=FALSE)
n_processed <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/D2/n_processed.tab", header=FALSE, stringsAsFactors=FALSE)
n_unique <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/D2/n_unique.tab", header=FALSE, stringsAsFactors=FALSE)

# quick visual comparison cortex vs liver
plot(n_unique$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Pseudomapping with Kallisto (D2 reference)", xlab="Samples", ylab="Uniquely mapped percentage",
     las=1)
legend("bottomright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")
plot(n_pseudoaligned$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Pseudomapping with Kallisto (D2 reference)", xlab="Samples", ylab="Pseudoaligned percentage",
     las=1)
legend("bottomright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")
plot(100-n_pseudoaligned$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Unmapped with Kallisto (D2 reference)", xlab="Samples", ylab="Unmapped percentage",
     las=1)
legend("topright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")


# stats with B6mm10 as a reference

# load stats
n_pseudoaligned <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_pseudoaligned.tab", header=FALSE, stringsAsFactors=FALSE)
n_processed <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_processed.tab", header=FALSE, stringsAsFactors=FALSE)
n_unique <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_unique.tab", header=FALSE, stringsAsFactors=FALSE)

# quick visual comparison cortex vs liver
plot(n_unique$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Pseudomapping with Kallisto (B6mm10 reference)", xlab="Samples", ylab="Uniquely mapped percentage",
     las=1)
legend("bottomright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")
plot(n_pseudoaligned$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Pseudomapping with Kallisto (B6mm10 reference)", xlab="Samples", ylab="Pseudoaligned percentage",
     las=1)
legend("bottomright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")
plot(100-n_pseudoaligned$V1/n_processed$V1*100, pch=19, col=as.factor(Meta$tissue),
     main="Unmapped with Kallisto (B6mm10 reference)", xlab="Samples", ylab="Unmapped percentage",
     las=1)
legend("topright", legend=levels(as.factor(Meta$tissue)), pch=19, col=1:2, title="Tissue")

#########################

# 16 December 2019

# Goal: plot D2 and B6 on the same plot

# load B6mm10 stats
n_pseudoaligned_B6mm10 <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_pseudoaligned.tab", header=FALSE, stringsAsFactors=FALSE)
n_processed_B6mm10 <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_processed.tab", header=FALSE, stringsAsFactors=FALSE)
n_unique_B6mm10 <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_unique.tab", header=FALSE, stringsAsFactors=FALSE)
# load D2 stats
n_pseudoaligned_D2 <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/D2/n_pseudoaligned.tab", header=FALSE, stringsAsFactors=FALSE)
n_processed_D2 <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/D2/n_processed.tab", header=FALSE, stringsAsFactors=FALSE)
n_unique_D2 <- read.table("F:/BXD/data/transcriptome/2.5_quantification_Kallisto/D2/n_unique.tab", header=FALSE, stringsAsFactors=FALSE)

# plot uniquely mapped
mylim <- c(min(n_unique_B6mm10$V1/n_processed_B6mm10$V1, n_unique_D2$V1/n_processed_D2$V1), max(n_unique_B6mm10$V1/n_processed_B6mm10$V1, n_unique_D2$V1/n_processed_D2$V1))*100
plot(n_unique_B6mm10$V1/n_processed_B6mm10$V1*100, col="black", bg="black",
     pch=19, ylim=mylim, las=1,
     main="Pseudomapping with Kallisto", xaxt="n", xlab="", ylab="Uniquely mapped percentage")
points(n_unique_D2$V1/n_processed_D2$V1*100, pch=19, col="tan", bg="tan")
for(myline in seq(0,2,0.5)){
        axis(1, at=86.5, line=myline, labels=c(""))
}
axis(1, at=c(43,129), line=0, tick=FALSE, labels=c("cortex", "liver"))
legend("topleft", legend=c("B6mm10","D2"), col=c("black","tan"), pch=19, title="Reference")

# plot pseudoaligned
mylim <- range(n_pseudoaligned_B6mm10$V1/n_processed_B6mm10$V1, n_pseudoaligned_D2$V1/n_processed_D2$V1)*100
plot(n_pseudoaligned_B6mm10$V1/n_processed_B6mm10$V1*100, col="black", bg="black",
     pch=19, ylim=mylim, las=1,
     main="Pseudomapping with Kallisto", xaxt="n", xlab="", ylab="Pseudoaligned percentage")
points(n_pseudoaligned_D2$V1/n_processed_D2$V1*100, pch=19, col="tan", bg="tan")
for(myline in seq(0,2,0.5)){
        axis(1, at=86.5, line=myline, labels=c(""))
}
axis(1, at=c(43,129), line=0, tick=FALSE, labels=c("cortex", "liver"))
legend("topleft", legend=c("B6mm10","D2"), col=c("black","tan"), pch=19, title="Reference")

# plot unmapped
mylim <- range(100-n_pseudoaligned_B6mm10$V1/n_processed_B6mm10$V1*100, 100-n_pseudoaligned_D2$V1/n_processed_D2$V1*100)
plot(100-n_pseudoaligned_B6mm10$V1/n_processed_B6mm10$V1*100, col="black", bg="black",
     pch=19, ylim=mylim, las=1,
     main="Unmapped with Kallisto", xaxt="n", xlab="", ylab="Unmapped percentage")
points(100-n_pseudoaligned_D2$V1/n_processed_D2$V1*100, pch=19, col="tan", bg="tan")
for(myline in seq(0,2,0.5)){
        axis(1, at=86.5, line=myline, labels=c(""))
}
axis(1, at=c(43,129), line=0, tick=FALSE, labels=c("cortex", "liver"))
legend("topleft", legend=c("B6mm10","D2"), col=c("black","tan"), pch=19, title="Reference")

# quick t-test for reference and for tissue
t.test(n_unique_B6mm10$V1/n_processed_B6mm10$V1*100, n_unique_D2$V1/n_processed_D2$V1*100, paired=TRUE)
t.test(n_unique_B6mm10$V1/n_processed_B6mm10$V1*100 ~ as.factor(Meta$tissue), paired=TRUE)
t.test(n_unique_D2$V1/n_processed_D2$V1*100 ~ as.factor(Meta$tissue), paired=TRUE)

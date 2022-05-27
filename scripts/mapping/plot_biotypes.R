# Biotype genes analysis of DM genes

# data
tot_genes <- c(4119, 2497, 13860, 11132)
protein_coding <- c(3970, 2442, 12967, 10794)

# plot
bp <- barplot(t(cbind(protein_coding,tot_genes-protein_coding)), pch=19, las=1, names.arg=c("DM_Cortex", "DM_Liver", "Total_Cortex", "Total_Liver"),col=c("grey40","grey70"))
text(bp, rep(1,4), label=paste(round(protein_coding/tot_genes*100),"%"), pos=3)

# GOAL: plot more biotypes
# load data
DM_cortexB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesDM_Cortex.tab", stringsAsFactors=FALSE)
DM_liverB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesDM_Liver.tab", stringsAsFactors=FALSE)
all_cortexB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesTotal_Cortex.tab", stringsAsFactors=FALSE)
all_liverB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesTotal_Liver.tab", stringsAsFactors=FALSE)
DM_cortexD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesDM_Cortex.tab", stringsAsFactors=FALSE)
DM_liverD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesDM_Liver.tab", stringsAsFactors=FALSE)
all_cortexD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesTotal_Cortex.tab", stringsAsFactors=FALSE)
all_liverD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesTotal_Liver.tab", stringsAsFactors=FALSE)

# adjusting margins
par(mar=c(4.1, 14.1, 2.1, 0.6))
# plotting
barplot(all_cortex$V1, names.arg=all_cortex$V2, las=1, horiz=TRUE, cex.names=0.9, xlab="Number of genes", main="Biotypes (Total Cortex)")
barplot(all_liver$V1, names.arg=all_liver$V2, las=1, horiz=TRUE, cex.names=0.9, xlab="Number of genes", main="Biotypes (Total Liver)")
barplot(DM_cortex$V1, names.arg=DM_cortex$V2, las=1, horiz=TRUE, cex.names=0.9, xlab="Number of genes", main="Biotypes (DM Cortex)")
barplot(DM_liver$V1, names.arg=DM_liver$V2, las=1, horiz=TRUE, cex.names=0.9, xlab="Number of genes", main="Biotypes (DM Liver)")

# Goal: plotting beside and use categories

# make a conversion table (biotypes to biotype categories)
biotypes <- unique(all_cortex$V2)
##c("protein_coding","pseudogenes","long_noncoding", "short_noncoding", "TEC")
categories <- c("long_noncoding", "long_noncoding", "protein_coding", "long_noncoding", "short_noncoding",
  "short_noncoding", "short_noncoding", "pseudogenes", "pseudogenes", "pseudogenes",      
  "protein_coding", "short_noncoding", "long_noncoding", "long_noncoding", "short_noncoding",
  "short_noncoding", "TEC", "pseudogenes", "pseudogenes", "pseudogenes", "protein_coding", 
  "pseudogenes")
ct <- cbind(biotypes,categories)

# convert
all_cortex$V3 <- ct[match(all_cortex$V2, ct[,1]),2]
all_liver$V3 <- ct[match(all_liver$V2, ct[,1]),2]
DM_cortex$V3 <- ct[match(DM_cortex$V2, ct[,1]),2]
DM_liver$V3 <- ct[match(DM_liver$V2, ct[,1]),2]

# sum genes by biotype categories
data <- array(NA, dim=c(5,4))
dimnames(data) <- list(unique(categories),c("DM_cortex", "DM_liver", "all_cortex", "all_liver"))

for(c in unique(categories)){
  data[c,"DM_cortex"] <- sum(DM_cortex$V1[DM_cortex$V3==c])
  data[c,"DM_liver"] <- sum(DM_liver$V1[DM_liver$V3==c])
  data[c,"all_cortex"] <- sum(all_cortex$V1[all_cortex$V3==c])
  data[c,"all_liver"] <- sum(all_liver$V1[all_liver$V3==c])
}

# adjusting margins
par(mar=c(3.1, 4.1, 3.1, 1.1))
barplot(data, las=1, beside=TRUE, col=gray.colors(5,0,1),
        ylab="Number of genes", main="Biotypes",
        legend.text=unique(categories), args.legend=list(x="topleft"))

barplot(t(t(data)/apply(data,2,sum))*100, ylim=c(0,100), las=1, beside=TRUE, col=gray.colors(5,0,1),
        ylab="% genes", main="Biotypes",
        legend.text=unique(categories), args.legend=list(x="topright"))

# 24 January 2020

# load data
DM_cortexB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesDM_Cortex.tab", stringsAsFactors=FALSE)
DM_liverB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesDM_Liver.tab", stringsAsFactors=FALSE)
all_cortexB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesTotal_Cortex.tab", stringsAsFactors=FALSE)
all_liverB6 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryB6BiotypesTotal_Liver.tab", stringsAsFactors=FALSE)
DM_cortexD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesDM_Cortex.tab", stringsAsFactors=FALSE)
DM_liverD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesDM_Liver.tab", stringsAsFactors=FALSE)
all_cortexD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesTotal_Cortex.tab", stringsAsFactors=FALSE)
all_liverD2 <- read.table("F:/BXD/data/transcriptome/6_Biotypes/summaryD2BiotypesTotal_Liver.tab", stringsAsFactors=FALSE)

# make a conversion table (biotypes to biotype categories)
biotypes <- unique(all_cortexB6$V2)
##c("protein_coding","pseudogenes","long_noncoding", "short_noncoding", "TEC")
categories <- c("long_noncoding", "long_noncoding", "protein_coding", "long_noncoding", "short_noncoding",
                "short_noncoding", "short_noncoding", "pseudogenes", "pseudogenes", "pseudogenes",      
                "protein_coding", "short_noncoding", "long_noncoding", "long_noncoding", "short_noncoding",
                "short_noncoding", "TEC", "pseudogenes", "pseudogenes", "pseudogenes", "protein_coding", 
                "pseudogenes")
ct <- cbind(biotypes,categories)

# convert
all_cortexB6$V3 <- ct[match(all_cortexB6$V2, ct[,1]),2]
all_liverB6$V3 <- ct[match(all_liverB6$V2, ct[,1]),2]
DM_cortexB6$V3 <- ct[match(DM_cortexB6$V2, ct[,1]),2]
DM_liverB6$V3 <- ct[match(DM_liverB6$V2, ct[,1]),2]
all_cortexD2$V3 <- ct[match(all_cortexD2$V2, ct[,1]),2]
all_liverD2$V3 <- ct[match(all_liverD2$V2, ct[,1]),2]
DM_cortexD2$V3 <- ct[match(DM_cortexD2$V2, ct[,1]),2]
DM_liverD2$V3 <- ct[match(DM_liverD2$V2, ct[,1]),2]

# sum genes by biotype categories
dataB6 <- array(NA, dim=c(5,4))
dimnames(dataB6) <- list(unique(categories),c("DM_cortex", "DM_liver", "all_cortex", "all_liver"))

for(c in unique(categories)){
  dataB6[c,"DM_cortex"] <- sum(DM_cortexB6$V1[DM_cortexB6$V3==c])
  dataB6[c,"DM_liver"] <- sum(DM_liverB6$V1[DM_liverB6$V3==c])
  dataB6[c,"all_cortex"] <- sum(all_cortexB6$V1[all_cortexB6$V3==c])
  dataB6[c,"all_liver"] <- sum(all_liverB6$V1[all_liverB6$V3==c])
}

dataD2 <- array(NA, dim=c(5,4))
dimnames(dataD2) <- list(unique(categories),c("DM_cortex", "DM_liver", "all_cortex", "all_liver"))

for(c in unique(categories)){
  dataD2[c,"DM_cortex"] <- sum(DM_cortexD2$V1[DM_cortexD2$V3==c])
  dataD2[c,"DM_liver"] <- sum(DM_liverD2$V1[DM_liverD2$V3==c])
  dataD2[c,"all_cortex"] <- sum(all_cortexD2$V1[all_cortexD2$V3==c])
  dataD2[c,"all_liver"] <- sum(all_liverD2$V1[all_liverD2$V3==c])
}

# adjusting margins
par(mar=c(3.1, 4.1, 3.1, 1.1))
# plotting numbers
barplot(dataB6, las=1, beside=TRUE, col=gray.colors(5,0,1),
        ylab="Number of genes", main="Biotypes (B6)",
        legend.text=unique(categories), args.legend=list(x="topleft"))
# plotting percentages
barplot(t(t(dataB6)/apply(dataB6,2,sum))*100, ylim=c(0,100), las=1, beside=TRUE, col=gray.colors(5,0,1),
        ylab="% genes", main="Biotypes (B6)",
        legend.text=unique(categories), args.legend=list(x="topright"))
# plotting numbers
barplot(dataD2, las=1, beside=TRUE, col=gray.colors(5,0,1),
        ylab="Number of genes", main="Biotypes (D2)",
        legend.text=unique(categories), args.legend=list(x="topleft"))
# plotting percentages
barplot(t(t(dataD2)/apply(dataD2,2,sum))*100, ylim=c(0,100), las=1, beside=TRUE, col=gray.colors(5,0,1),
        ylab="% genes", main="Biotypes (D2)",
        legend.text=unique(categories), args.legend=list(x="topright"))

# 27 May 2022

# calculate enrichment
round(dataB6["pseudogenes","DM_cortex"] / sum(dataB6[,"DM_cortex"]) / (dataB6["pseudogenes","all_cortex"] / sum(dataB6[,"all_cortex"])), 1)
round(dataB6["pseudogenes","DM_liver"] / sum(dataB6[,"DM_liver"]) / (dataB6["pseudogenes","all_liver"] / sum(dataB6[,"all_liver"])), 1)
round(dataD2["pseudogenes","DM_cortex"] / sum(dataD2[,"DM_cortex"]) / (dataD2["pseudogenes","all_cortex"] / sum(dataD2[,"all_cortex"])), 1)
round(dataD2["pseudogenes","DM_liver"] / sum(dataD2[,"DM_liver"]) / (dataD2["pseudogenes","all_liver"] / sum(dataD2[,"all_liver"])), 1)

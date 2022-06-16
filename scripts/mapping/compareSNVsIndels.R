# compare for D2 samples reference with SNVs vs indels vs both

# set working directory
setwd("F:/BXD/analysis/Mapping")

#load needed libraries
library(knitr)
library(vioplot)

# source function
source("F:/BXD/analysis/scripts/mapping/retrieveSTARstatistics.R")

path <- "F:/BXD/data/transcriptome/2_Mapping_STAR/OnPersonalizedRef/"
statsName <- "Uniquely mapped reads % |"
references <- dir(path, pattern="DBA_2J_genome_")
ref <- "DBA_2J_genome_nonrandomized_indels_maternal"
samples <- gsub("_Log.final.out","", list.files(paste0(path, "/", ref), pattern="_Log.final.out"))

resUniqP <- outer(FUN=retrieveStats, X=references, Y=paste0(samples, "_"), path="F:/BXD/data/transcriptome/2_Mapping_STAR/OnPersonalizedRef/", statsName="Uniquely mapped reads % |")
rownames(resUniqP) <- references
colnames(resUniqP) <- samples

baseline <- retrieveStats("B6mm10_primaryassembly", paste0(samples, "_exactunique"), "F:/BXD/data/transcriptome/2_Mapping_STAR/OnGenome", statsName)
resUniqP <- rbind(resUniqP, baseline)

# selection of references and samples of interest
row_select <- intersect(grep("_randomized_", rownames(resUniqP), value=TRUE, invert=TRUE),grep("maternal", rownames(resUniqP), value=TRUE, invert=TRUE))
col_select <- grep("(DB)(1|2)", colnames(resUniqP), value=TRUE)
data <- resUniqP[row_select, col_select]

# formatting (simplification of labels)
rownames(data) <- gsub(pattern="DBA_2J_genome_nonrandomized_", replacement="", rownames(data))
rownames(data) <- gsub(pattern="_paternal", replacement="", rownames(data))
rownames(data) <- gsub(pattern="baseline", replacement="GRCm38 only", rownames(data))

# plot each sample
par(mfrow=c(2,2), mar=c(7.1, 4.1, 2.1, 2.1))
samples <- colnames(data)
for(sample in samples){
  barplot(sort(data[(-c(1,5)), sample]) - data["GRCm38 only", sample], col="tan",
          main=sample, las=2, cex.names=0.9, ylab="Difference in Uniquely mapped reads %", ylim=c(-0,7))
}

# one plot for all samples
par(mfrow=c(1,1))
plot(0:2, 0:2, ylim=c(0,6.5), col="white", las=1, xaxt="n", xlab="Reference modified with D2-specific variants", ylab="Difference in Uniquely mapped reads %")
axis(1, at=0:2, labels=c("indels", "SNVs", "indelsSNVs"))
points(rep(0,4), data["indels", ] - data["GRCm38 only", ], pch=21)
points(rep(1,4), data["SNVs", ] - data["GRCm38 only", ], pch=21)
points(rep(2,4), data["indelsSNVs", ] - data["GRCm38 only", ], pch=21)

# average gain per variant type
round(mean(data["indels", ] - data["GRCm38 only", ]),1)
round(mean(data["SNVs", ] - data["GRCm38 only", ]),1)
round(mean(data["indelsSNVs", ] - data["GRCm38 only", ]),1)

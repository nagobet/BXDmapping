# Retrieve gene position

# GOAL: retrieve gene position from gtf file

# DATE: 20200624

# adapted from Get_Gene_Position.Rmd by Maxime

# INPUT
# 1 file with gene annotation (.gtf)
# OUTPUT
# 1 file with gene name and positions (.txt)

#############################################################

# load gtf data
gtf <- read.table("/mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf", stringsAsFactors=FALSE, sep="\t")
colnames(gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Frame", "info")

# formate data
# extract only features of type gene (filter out transcripts, exons, introns)
genegtf <- subset(gtf, subset=(gtf$Type=="gene"))
# add "chr" to chromosome
genegtf$Chromosome <- paste0("chr", genegtf$Chromosome)
# isolate gene name
genegtf$GeneName <- gsub(" gene_name ", "", grep("gene_name", unlist(strsplit(genegtf$info, ";")), value=TRUE))
# calculate mean, start, and end position of gene in Mb
genegtf$Position <- round(apply(genegtf[,c("Start","End")], MARGIN=1, FUN=mean)/1e6, 6)
genegtf$Start <- round(genegtf$Start/1e6, 6)
genegtf$End <- round(genegtf$End/1e6, 6)
# keep columns of interest and remove duplicated gene names
genepos <- subset(genegtf, subset=!duplicated(genegtf$GeneName), select=c("Chromosome", "Position", "Start", "End"))
rownames(genepos) <- genegtf$GeneName[!duplicated(genegtf$GeneName)]

# save output
write.table(genepos, file="/mnt/nas/BXD/references/transcriptome/GenePosition.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

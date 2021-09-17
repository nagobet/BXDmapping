#!/usr/bin/env Rscript
# Gene counts normalization

# GOAL: normalize gene counts (tpm)

# DATE: 20200616 to 20200619

# install packages if necessary
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if(!require("limma",character.only=TRUE)) BiocManager::install("limma")
if(!require("edgeR",character.only=TRUE)) BiocManager::install("edgeR")

# load needed libraries
library(limma)
library(edgeR)

# load table to convert mouse line names
##LinesNames <- read.table("F:/BXD/data/ConvertLineNames.tsv", header=FALSE, stringsAsFactors=FALSE, sep="\t")
LinesNames <- read.table("/mnt/nas/BXD/data/ConvertLineNames.tsv", header=FALSE, stringsAsFactors=FALSE, sep="\t")
# load table to convert gene id to gene name
##GeneConvert <- read.table("F:/BXD/references/gene_names_B6mm10.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
GeneConvert <- read.table("/mnt/nas/BXD/references/gene_names_B6mm10.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")

# INPUT
# one file with all counts (condition and tissues together).
# OUTPUT
# 8 files of normalized tpm counts (separated by tissue and condition), log2 or not.

#####################################

# handle Rscript arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Error: argument missing/incorrect (please provide path to Summary sample count file, ending with /)", call.=FALSE)
}
path <- args[1]
print(path)
print("tpm normalization")

# load counts data, formate, and split by tissue
all_counts <- read.table(paste0(path, "Summary_ReadsPerGene.out.tab"), sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names="GeneID")
colnames(all_counts) <- gsub("X", "C", colnames(all_counts))
idx_gene <- grep("ENSMUSG", rownames(all_counts), value=FALSE)
idx_uniq <- which(!duplicated(GeneConvert$V2[match(rownames(all_counts), GeneConvert$V1)]))
idx_keep <- intersect(idx_gene, idx_uniq)
all_counts <- all_counts[idx_keep,]
rownames(all_counts) <- GeneConvert$V2[match(rownames(all_counts), GeneConvert$V1)]
cortex_counts <- subset(all_counts, select=grep("C", colnames(all_counts)))
liver_counts <- subset(all_counts, select=grep("L", colnames(all_counts)))

# Normalize by gene length
# prepare gene length for different samples
linenumberbysample <- as.numeric(gsub("[CLnsd]", "", colnames(cortex_counts)))
gene_length_matrix <- matrix(NA, nrow=nrow(cortex_counts), ncol=ncol(cortex_counts))
# load gene length
##gl <- read.table(paste0("F:/BXD/references/transcriptome/genelength.txt"), header=TRUE, stringsAsFactors=FALSE, row.names="gene")
gl <- read.table(paste0("/mnt/nas/BXD/references/transcriptome/genelength.txt"), header=TRUE, stringsAsFactors=FALSE, row.names="gene")
p <- 1
for(l in linenumberbysample){
    # retrieve merged gene length in sample genes order
  gene_length_matrix[,p] <- gl[GeneConvert$V1[match(rownames(cortex_counts), GeneConvert$V2)],"merged"]
  p <- p + 1
}
# normalize by gene length in kilobases
cortex_tpm <- cortex_counts/(gene_length_matrix/1000)
liver_tpm <- liver_counts/(gene_length_matrix/1000)

# filter lowly expressed genes
dC <- DGEList(counts=cortex_tpm)
dL <- DGEList(counts=liver_tpm)

filt_logical <- rowSums(cpm(dC)>0.5)>=20
dCf <- dC[filt_logical,]
rownames(dCf) <- rownames(dC)[filt_logical]
dim(dCf)

filt_logical <- rowSums(cpm(dL)>0.5)>=20
dLf <- dL[filt_logical,]
rownames(dLf) <- rownames(dL)[filt_logical]
dim(dLf)

# compute normalization factors
dCf <- calcNormFactors(dCf, method="TMM")
dCf <- estimateCommonDisp(dCf, verbose=TRUE)
dLf <- calcNormFactors(dLf, method="TMM")
dLf <- estimateCommonDisp(dLf, verbose=TRUE)

# split by condition (nsd or sd)
cnsd <- dCf[, grep("nsd", colnames(dCf))]
csd <- dCf[, grep("nsd", colnames(dCf), invert=TRUE)]
lnsd <- dLf[, grep("nsd", colnames(dLf))]
lsd <- dLf[, grep("nsd", colnames(dLf), invert=TRUE)]

# compute TMM normalized TPM counts
cnsdtpm <- data.frame(cpm(cnsd, normalized.lib.sizes=TRUE))
csdtpm <- data.frame(cpm(csd, normalized.lib.sizes=TRUE))
lnsdtpm <- data.frame(cpm(lnsd, normalized.lib.sizes=TRUE))
lsdtpm <- data.frame(cpm(lsd, normalized.lib.sizes=TRUE))

# reformate column names with BXD lines
convertColNames <- function(data){
  linenumber <- gsub("[CLnsd]", "", colnames(data))
  return(LinesNames$V5[match(linenumber, LinesNames$V10)])
}
colnames(cnsdtpm) <- convertColNames(cnsdtpm)
colnames(csdtpm) <- convertColNames(csdtpm)
colnames(lnsdtpm) <- convertColNames(lnsdtpm)
colnames(lsdtpm) <- convertColNames(lsdtpm)

# save TMM normalized TPM counts into files
write.table(cnsdtpm, file=paste0(path, "TMMnormalized_TPM_Cortex_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(csdtpm, file=paste0(path, "TMMnormalized_TPM_Cortex_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lnsdtpm, file=paste0(path, "TMMnormalized_TPM_Liver_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lsdtpm, file=paste0(path, "TMMnormalized_TPM_Liver_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)

# compute log2 normalized counts
cnsdtpmlog <- data.frame(cpm(cnsd, TMMnormalized.lib.sizes=TRUE, log=TRUE, prior.count=1))
csdtpmlog <- data.frame(cpm(csd, TMMnormalized.lib.sizes=TRUE, log=TRUE, prior.count=1))
lnsdtpmlog <- data.frame(cpm(lnsd, TMMnormalized.lib.sizes=TRUE, log=TRUE, prior.count=1))
lsdtpmlog <- data.frame(cpm(lsd, TMMnormalized.lib.sizes=TRUE, log=TRUE, prior.count=1))

# reformate column names with BXD lines
colnames(cnsdtpmlog) <- convertColNames(cnsdtpmlog)
colnames(csdtpmlog) <- convertColNames(csdtpmlog)
colnames(lnsdtpmlog) <- convertColNames(lnsdtpmlog)
colnames(lsdtpmlog) <- convertColNames(lsdtpmlog)

# save log2 normalized counts into files
write.table(cnsdtpmlog, file=paste0(path, "TMMnormalized_log2TPM_Cortex_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(csdtpmlog, file=paste0(path, "TMMnormalized_log2TPM_Cortex_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lnsdtpmlog, file=paste0(path, "TMMnormalized_log2TPM_Liver_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lsdtpmlog, file=paste0(path, "TMMnormalized_log2TPM_Liver_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)

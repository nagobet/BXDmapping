#!/usr/bin/env Rscript
# Gene counts normalization

# GOAL: normalize gene counts (cpm)

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
# 8 files of normalized cpm counts (separated by tissue and condition), log2 or not.

#####################################

# handle Rscript arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Error: argument missing/incorrect (please provide path to Summary sample count file, ending with /)", call.=FALSE)
}
path <- args[1]
##path <- "F:/BXD/data/transcriptome/2_mapping_STAR/MappingParametersOptimization/genotypesandimputed_withannotation_EndToEnd_0_0/"
print(path)
print("cpm normalization")

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

# filter lowly expressed genes
dC <- DGEList(counts=cortex_counts)
dL <- DGEList(counts=liver_counts)

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

# compute TMM normalized CPM counts
cnsdcpm <- data.frame(cpm(cnsd, normalized.lib.sizes=TRUE))
csdcpm <- data.frame(cpm(csd, normalized.lib.sizes=TRUE))
lnsdcpm <- data.frame(cpm(lnsd, normalized.lib.sizes=TRUE))
lsdcpm <- data.frame(cpm(lsd, normalized.lib.sizes=TRUE))

# reformate column names with BXD lines
convertColNames <- function(data){
  linenumber <- gsub("[CLnsd]", "", colnames(data))
  return(LinesNames$V5[match(linenumber, LinesNames$V10)])
}
colnames(cnsdcpm) <- convertColNames(cnsdcpm)
colnames(csdcpm) <- convertColNames(csdcpm)
colnames(lnsdcpm) <- convertColNames(lnsdcpm)
colnames(lsdcpm) <- convertColNames(lsdcpm)

# save TMM normalized CPM counts into files
write.table(cnsdcpm, file=paste0(path, "TMMnormalized_CPM_Cortex_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(csdcpm, file=paste0(path, "TMMnormalized_CPM_Cortex_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lnsdcpm, file=paste0(path, "TMMnormalized_CPM_Liver_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lsdcpm, file=paste0(path, "TMMnormalized_CPM_Liver_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)

# compute log2 normalized counts
cnsdcpmlog <- data.frame(cpm(cnsd, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1))
csdcpmlog <- data.frame(cpm(csd, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1))
lnsdcpmlog <- data.frame(cpm(lnsd, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1))
lsdcpmlog <- data.frame(cpm(lsd, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1))

# reformate column names with BXD lines
colnames(cnsdcpmlog) <- convertColNames(cnsdcpmlog)
colnames(csdcpmlog) <- convertColNames(csdcpmlog)
colnames(lnsdcpmlog) <- convertColNames(lnsdcpmlog)
colnames(lsdcpmlog) <- convertColNames(lsdcpmlog)

# save log2 normalized counts into files
write.table(cnsdcpmlog, file=paste0(path, "TMMnormalized_log2CPM_Cortex_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(csdcpmlog, file=paste0(path, "TMMnormalized_log2CPM_Cortex_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lnsdcpmlog, file=paste0(path, "TMMnormalized_log2CPM_Liver_NSD.tab"), row.names=TRUE, sep="\t", quote=FALSE)
write.table(lsdcpmlog, file=paste0(path, "TMMnormalized_log2CPM_Liver_SD.tab"), row.names=TRUE, sep="\t", quote=FALSE)

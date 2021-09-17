# Testing filtering

# GOAL: Perform different filtering of lowly expressed genes and normalization

# DATE: 30 July 2020 - 31 August 2020

# load needed libraries
library(limma)
library(edgeR)
library(zoo)
library(sm)
library(vioplot)

# Variables
inpath <- "F:/BXD/data/transcriptome/2_mapping_STAR/MappingParametersOptimization/genotypesandimputed_withannotation_Local_0_2/"
outpath <- "F:/BXD/data/transcriptome/4_normalization/FilteringTesting/"
LinesNames <- read.table("F:/BXD/data/ConvertLineNames.tsv", header=FALSE, stringsAsFactors=FALSE, sep="\t")
GeneConvert <- read.table("F:/BXD/references/gene_names_B6mm10.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")

###############

# function for filtering and normalization (thresholds)
filternormalizeCounts_thresholds <- function(min_cpm, min_samples){
  
  # load counts data, formate, and split by tissue
  all_counts <- read.table(paste0(inpath, "Summary_ReadsPerGene.out.tab"), sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names="GeneID")
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
  
  filt_logical <- rowSums(cpm(dC)>min_cpm)>=min_samples
  dCf <- dC[filt_logical,]
  rownames(dCf) <- rownames(dC)[filt_logical]
  dim(dCf)
  
  filt_logical <- rowSums(cpm(dL)>min_cpm)>=min_samples
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
  write.table(cnsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Cortex_NSD_", min_cpm, "_", min_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(csdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Cortex_SD_", min_cpm, "_", min_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(lnsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Liver_NSD_", min_cpm, "_", min_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(lsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Liver_SD_", min_cpm, "_", min_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)

  min_cpm <- gsub(pattern=".", replacement="p", min_cpm, fixed=TRUE)
  print(paste("Filtering and normalization for", min_cpm, "cpm and", min_samples, "samples is DONE!"))
}

# sufficient samples with sufficient cpm (thresholds)
filternormalizeCounts_thresholds(min_cpm=0.5, min_samples=20)
filternormalizeCounts_thresholds(min_cpm=0.5, min_samples=44)
filternormalizeCounts_thresholds(min_cpm=0.5, min_samples=22)
filternormalizeCounts_thresholds(min_cpm=0.5, min_samples=11)
filternormalizeCounts_thresholds(min_cpm=1, min_samples=20)
filternormalizeCounts_thresholds(min_cpm=0.1, min_samples=20)

###############

# function for filtering and normalization (thresholds)
filternormalizeCounts_average <- function(percentile){
  
  # load counts data, formate, and split by tissue
  all_counts <- read.table(paste0(inpath, "Summary_ReadsPerGene.out.tab"), sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names="GeneID")
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
  
  threshold <- quantile(apply(dC$counts, MARGIN=1, FUN=median), probs=percentile/100)
  filt_logical <- apply(dC$counts, MARGIN=1, FUN=median)>=threshold
  dCf <- dC[filt_logical,]
  rownames(dCf) <- rownames(dC)[filt_logical]
  dim(dCf)
  
  threshold <- quantile(apply(dL$counts, MARGIN=1, FUN=median), probs=percentile/100)
  filt_logical <- apply(dL$counts, MARGIN=1, FUN=median)>=threshold
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
  write.table(cnsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Cortex_NSD_", percentile, ".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(csdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Cortex_SD_", percentile, ".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(lnsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Liver_NSD_",percentile, ".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(lsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Liver_SD_", percentile, ".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  
  print(paste("Filtering and normalization for", percentile, "% is DONE!"))
}

# minimal average gene counts (average)
filternormalizeCounts_average(percentile=0)
filternormalizeCounts_average(percentile=75)

###############

# create a function for filtering and normalization (not0) only genes where at least number_samples samples have a non 0 gene count
filternormalizeCounts_not0 <- function(number_samples){
  
  # load counts data, formate, and split by tissue
  all_counts <- read.table(paste0(inpath, "Summary_ReadsPerGene.out.tab"), sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names="GeneID")
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
  
  filt_logical <- rowSums(dC$counts>0)>=number_samples
  dCf <- dC[filt_logical,]
  rownames(dCf) <- rownames(dC)[filt_logical]
  dim(dCf)
  
  filt_logical <- rowSums(dC$counts>0)>=number_samples
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
  write.table(cnsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Cortex_NSD_", "not0", "_", number_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(csdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Cortex_SD_", "not0", "_", number_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(lnsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Liver_NSD_", "not0", "_", number_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  write.table(lsdcpmlog, file=paste0(outpath, "TMMnormalized_log2CPM_Liver_SD_", "not0", "_", number_samples,".tab"), row.names=TRUE, sep="\t", quote=FALSE)
  
  print(paste("Filtering and normalization for at least", number_samples, "of samples with more than 0 count is DONE!"))
}

# Use different limit number of samples
filternormalizeCounts_not0(33)
filternormalizeCounts_not0(16)
filternormalizeCounts_not0(49)
filternormalizeCounts_not0(66)

###############

# test if normalized counts previously computed (on Linux, "old") are the same that the ones just computed (on Windows, "new") with the same filtering
old <- read.table(paste0(inpath, "TMMnormalized_log2CPM_Cortex_NSD.tab"), header=TRUE, stringsAsFactors=FALSE)
new <- read.table(paste0(outpath, "TMMnormalized_log2CPM_Cortex_NSD_0.5_20.tab"), header=TRUE, stringsAsFactors=FALSE)
for(col in colnames(new)[-1]){
  plot(old[,col]-new[,col], main=col)
}
#It is not, but difference is very low (maximum 1e-13).

# get session details
sessionInfo()

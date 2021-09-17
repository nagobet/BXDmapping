# TMM normalization diagnostic plot

# Goal: check effect of TMM normalization

# source: https://www.biostars.org/p/9465851/#9465854

# load needed libraries
library(limma)
library(edgeR)

# load table to convert mouse line names
LinesNames <- read.table("F:/BXD/data/ConvertLineNames.tsv", header=FALSE, stringsAsFactors=FALSE, sep="\t")
# load table to convert gene id to gene name
GeneConvert <- read.table("F:/BXD/references/gene_names_B6mm10.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")

# load counts data, formate, and split by tissue
all_counts <- read.table("F:/BXD/data/MappingEvaluation/genotypesandimputed_withannotation_Local_0_10/Summary_ReadsPerGene.out.tab", sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names="GeneID")
colnames(all_counts) <- gsub("X", "C", colnames(all_counts))
idx_gene <- grep("ENSMUSG", rownames(all_counts), value=FALSE)
idx_uniq <- which(!duplicated(GeneConvert$V2[match(rownames(all_counts), GeneConvert$V1)]))
idx_keep <- intersect(idx_gene, idx_uniq)
all_counts <- all_counts[idx_keep,]
rownames(all_counts) <- GeneConvert$V2[match(rownames(all_counts), GeneConvert$V1)]
cortex_counts <- subset(all_counts, select=grep("C", colnames(all_counts)))
liver_counts <- subset(all_counts, select=grep("L", colnames(all_counts)))

dC <- DGEList(counts=cortex_counts)
dL <- DGEList(counts=liver_counts)

# compute unfiltered not normalized log2 CPM counts
ccpm_unfiltered <- edgeR::cpmByGroup(dC, group=grepl("nsd", colnames(dC)), log=TRUE)
lcpm_unfiltered <- edgeR::cpmByGroup(dL, group=grepl("nsd", colnames(dL)), log=TRUE)

# filter lowly expressed genes
filt_logical <- rowSums(cpm(dC)>0.5)>=20
dCf <- dC[filt_logical,]
rownames(dCf) <- rownames(dC)[filt_logical]
dim(dCf)

filt_logical <- rowSums(cpm(dL)>0.5)>=20
dLf <- dL[filt_logical,]
rownames(dLf) <- rownames(dL)[filt_logical]
dim(dLf)

# compute not normalized log2 CPM counts
ccpm <- edgeR::cpmByGroup(dCf, group=grepl("nsd", colnames(dCf)), log=TRUE)
lcpm <- edgeR::cpmByGroup(dLf, group=grepl("nsd", colnames(dLf)), log=TRUE)

# compute normalization factors
dCf <- calcNormFactors(dCf, method="TMM")
##dCf <- estimateCommonDisp(dCf, verbose=TRUE)
dLf <- calcNormFactors(dLf, method="TMM")
##dLf <- estimateCommonDisp(dLf, verbose=TRUE)

# compute TMM normalized log2 CPM counts
ccpm_TMM <- edgeR::cpmByGroup(dCf, group=grepl("nsd", colnames(dCf)), log=TRUE)
lcpm_TMM <- edgeR::cpmByGroup(dLf, group=grepl("nsd", colnames(dLf)), log=TRUE)

# Make the MA-plot using smoothScatter:
par(mfrow=c(1,3))

myxlim=c(0,500)
myylim=c(-10,10)

# no filtering cortex
g1 <- ccpm_unfiltered[,1]
g2 <- ccpm_unfiltered[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "unfiltered in cortex",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

# naive cortex
g1 <- ccpm[,1]
g2 <- ccpm[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "naive in cortex",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

# TMM cortex
g1 <- ccpm_TMM[,1]
g2 <- ccpm_TMM[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "TMM in cortex",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

# no filtering liver
g1 <- lcpm_unfiltered[,1]
g2 <- lcpm_unfiltered[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "unfiltered in liver",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

# naive liver
g1 <- lcpm[,1]
g2 <- lcpm[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "naive in liver",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

# TMM liver
g1 <- lcpm_TMM[,1]
g2 <- lcpm_TMM[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "TMM in liver",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

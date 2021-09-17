# TMM normalization diagnostic plot

# Goal: check effect of TMM normalization on tissues (normalize the 2 tissues together)

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

# filter lowly expressed genes
dall <- DGEList(counts=all_counts)

filt_logical <- rowSums(cpm(dall)>0.5)>=20
dallf <- dall[filt_logical,]
rownames(dallf) <- rownames(dall)[filt_logical]
dim(dallf)

# compute not normalized log2 CPM counts
allcpm <- edgeR::cpmByGroup(dallf, group=grepl("C", colnames(dallf)), log=TRUE)

# compute normalization factors
dalldCf <- calcNormFactors(dallf, method="TMM")

# compute TMM normalized log2 CPM counts
allcpm_TMM <- edgeR::cpmByGroup(dalldCf, group=grepl("C", colnames(dalldCf)), log=TRUE)

# Make the MA-plot using smoothScatter:
par(mfrow=c(1,2))

myxlim=c(0,100)
myylim=c(-10,10)

# naive in both tissues
g1 <- allcpm[,1]
g2 <- allcpm[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "naive in both tissues",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

# TMM in both tissues
g1 <- allcpm_TMM[,1]
g2 <- allcpm_TMM[,2]
smoothScatter(x = 0.5 * (g1+g2), # average expr
              y = g1-g2, # fold change
              xlab = "average log expression",
              ylab = "log fold change",
              main = "TMM in both tissues",
              xlim=myxlim,
              ylim=myylim)
abline(h=0, col="red")
abline(h=log2(2), col="red", lty=2)
abline(h=-log2(2), col="red", lty=2)

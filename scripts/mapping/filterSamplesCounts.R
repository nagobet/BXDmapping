# filter Samples Counts

# GOAL: filter out old BXD (05, 29, 29t, 32), problematic BXD (63), parental, and F1 samples counts

# DATE: 20200626

# INPUT: 1 summary file of gene counts for all samples
# OUTPUT: 1 summary file of gene counts for BXD samples

counts <- read.table("/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/SummaryAll_ReadsPerGene.out.tab", stringsAsFactors=FALSE, header=TRUE)
##counts <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/SummaryAll_ReadsPerGene.out.tab", stringsAsFactors=FALSE, header=TRUE)

keep_names <- grep("05|29|32|63|B|D.", colnames(counts), invert=TRUE, value=TRUE)
countsFiltered <- subset(counts, select=keep_names)
  
write.table(countsFiltered, file="/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/Summary_ReadsPerGene.out.tab", quote=FALSE, sep="\t", row.names=FALSE)
##write.table(countsFiltered, file="F:/BXD/data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/Summary_ReadsPerGene.out.tab", quote=FALSE, sep="\t", row.names=FALSE)

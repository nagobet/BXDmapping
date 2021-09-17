# filter Samples Counts

# GOAL: filter out old BXD (05, 29, 29t, 32), problematic BXD (63), parental, and F1 samples counts

# DATE: 20201208

# INPUT: 1 summary file of gene counts for all samples
# OUTPUT: 1 summary file of gene counts for BXD samples

counts <- read.table("/mnt/nas/BXD/data/MappingEvaluation/genotypesandimputed_withoutannotation_EndToEnd_1_0/SummaryAll_ReadsPerGene.out.tab", stringsAsFactors=FALSE, header=TRUE)
##counts <- read.table("F:/BXD/data/MappingEvaluation/genotypesandimputed_withoutannotation_EndToEnd_1_0/SummaryAll_ReadsPerGene.out.tab", stringsAsFactors=FALSE, header=TRUE)

keep_names <- grep("05|29|32|63|B|D.", colnames(counts), invert=TRUE, value=TRUE)
countsFiltered <- subset(counts, select=keep_names)
  
write.table(countsFiltered, file="/mnt/nas/BXD/data/MappingEvaluation/genotypesandimputed_withoutannotation_EndToEnd_1_0/Summary_ReadsPerGene.out.tab", quote=FALSE, sep="\t", row.names=FALSE)
##write.table(countsFiltered, file="F:/BXD/data/MappingEvaluation/genotypesandimputed_withoutannotation_EndToEnd_1_0/Summary_ReadsPerGene.out.tab", quote=FALSE, sep="\t", row.names=FALSE)

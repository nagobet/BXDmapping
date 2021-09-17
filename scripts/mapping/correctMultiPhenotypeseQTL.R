# cis-eQTL multi test (phenotypes) correction

# GOAL: compute q-value for eQTL to assess that multiple phenotypes (genes expression) were tested.

# DATE: 20200623

# INPUT
# 1 file with list of eQTLs from FastQTL
# OUTPUT
# 1 file with list of q-value corrected eQTLs

#####################################

# load needed library
library(qvalue)

# handle Rscript arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Error: argument missing/incorrect (please file name with its path, but not .txt)", call.=FALSE)
}
eQTLfile <- args[1]
##path <- "F:/BXD/data/transcriptome/2_mapping_STAR/MappingParametersOptimization/genotypesandimputed_withannotation_EndToEnd_0_0/"
print(eQTLfile)
print("cis-eQTL multi test (phenotypes) correction")

eQTLdata <- read.table(paste0(eQTLfile, ".txt"), stringsAsFactors=FALSE)

# remove a duplicated gene Cdkn2d
eQTLdata <- eQTLdata[!duplicated(eQTLdata$V1),]

# selects column with p-value (adjusted for different variants tested for this gene) and id of the best variant
eQTLdataqvalue <- eQTLdata[,c("V11","V6")]
eQTLdataqvalue$V11 <- qvalue(eQTLdata$V11)$qvalues
rownames(eQTLdataqvalue) <- eQTLdata$V1
colnames(eQTLdataqvalue) <- c("adjustedpvalue","marker")
eQTLdataqvalue$adjustedpvalue[is.na(eQTLdataqvalue$adjustedpvalue)] <- 1

# save corrected p-values (q-values)
write.table(eQTLdataqvalue, file=paste(eQTLfile, "pvalcorrected.txt", sep="_"), quote=FALSE, row.names=TRUE, col.names=T, sep="\t")

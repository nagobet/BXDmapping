# set working directory
setwd("F:/BXD/analyses/Comparison_mapping_B6vsD2")

# load data
sam541_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151541_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam541_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151541_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)

sam543_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151543_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam543_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151543_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)

sam547_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151547_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam547_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151547_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)

sam549_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151549_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam549_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151549_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)


sam627_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151627_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam627_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151627_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)

sam629_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151629_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam629_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151629_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)

sam633_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151633_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam633_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151633_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)

sam635_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151635_uniqtoB6.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)
sam635_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/GSM3151635_uniqtoD2.samshort", comment.char="", header=FALSE, stringsAsFactors=FALSE)


##sam541_B6 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/tmp", comment.char="", header=FALSE, stringsAsFactors=FALSE)
##sam541_D2 <- read.table("../../data/transcriptome/2_mapping_STAR/D2vsB6/tmpD2", comment.char="", header=FALSE, stringsAsFactors=FALSE)


head(sam541_B6)

summary(sam541_B6)

summary(as.factor(sam541_B6$V2))
summary(as.factor(sam541_D2$V2))

summary(as.factor(sam543_B6$V2))
summary(as.factor(sam543_D2$V2))

summary(as.factor(sam547_B6$V2))
summary(as.factor(sam547_D2$V2))

summary(as.factor(sam549_B6$V2))
summary(as.factor(sam549_D2$V2))


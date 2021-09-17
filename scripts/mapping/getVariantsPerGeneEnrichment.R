# variants per gene enrichment in DM genes

# set working directory
setwd("F:/BXD/data/transcriptome/5_DifferentialMapping")

# checking for liver BXD vs B6

# load genes-variants overlap
lgS <- read.table("liver_genesSNVsoverlap_BXDvsB6mm10permissive.txt", header=FALSE, stringsAsFactors=FALSE)
lgi <- read.table("liver_genesindelsoverlap_BXDvsB6mm10permissive.txt", header=FALSE, stringsAsFactors=FALSE)

# load DM data
DM_l <- read.table("DifferentialMappingBXDvsB6mm10permissive_Limma_Liver.txt", header=TRUE, stringsAsFactors=FALSE)

# retrieve gene names of DM and not DM
genenames_sig <- rownames(DM_l)[which(DM_l$adj.P.Val<0.05)]
genenames_notsig <- rownames(DM_l)[which(DM_l$adj.P.Val>=0.05)]

# check if samples number is balanced
length(genenames_sig)
length(genenames_notsig)
#no

# normalize by the length of the genes
norm_sig <- lgS$V24[match(genenames_sig, lgS$V16)]/(lgS$V5[match(genenames_sig, lgS$V16)]-lgS$V4[match(genenames_sig, lgS$V16)])
norm_notsig <- lgS$V24[match(genenames_notsig, lgS$V16)]/(lgS$V5[match(genenames_notsig, lgS$V16)]-lgS$V4[match(genenames_notsig, lgS$V16)])

summary(norm_sig)
summary(norm_notsig)
sd(norm_sig)
sd(norm_notsig)

# check if normally distributed
hist(norm_sig, breaks=50)
hist(norm_notsig, breaks=50)
#no, not even symetric.

t.test(norm_sig, norm_notsig)
# Average 48 SNVs by kb of DM genes, average 18 SNVs by kb of non-DM genes. 

#################################################################3

# generalizing

# create a function check if there is an enrichment
testEnrichment <- function(DM_data, genevariantsoverlap){
  # retrieve gene names of DM and not DM
  genenames_sig <- rownames(DM_data)[which(DM_data$adj.P.Val<0.05)]
  genenames_notsig <- rownames(DM_data)[which(DM_data$adj.P.Val>=0.05)]
  
  # normalize by the length of the genes
  norm_sig <- genevariantsoverlap$V24[match(genenames_sig, genevariantsoverlap$V16)]/(genevariantsoverlap$V5[match(genenames_sig, genevariantsoverlap$V16)]-genevariantsoverlap$V4[match(genenames_sig, genevariantsoverlap$V16)])
  norm_notsig <- genevariantsoverlap$V24[match(genenames_notsig, genevariantsoverlap$V16)]/(genevariantsoverlap$V5[match(genenames_notsig, genevariantsoverlap$V16)]-genevariantsoverlap$V4[match(genenames_notsig, genevariantsoverlap$V16)])
  
  t.test(norm_sig, norm_notsig)
}

# load genes-variants overlap
lgS <- read.table("liver_genesSNVsoverlap_BXDvsB6mm10permissive.txt", header=FALSE, stringsAsFactors=FALSE)
lgi <- read.table("liver_genesindelsoverlap_BXDvsB6mm10permissive.txt", header=FALSE, stringsAsFactors=FALSE)

# load DM data
DM_l <- read.table("DifferentialMappingBXDvsB6mm10permissive_Limma_Liver.txt", header=TRUE, stringsAsFactors=FALSE)

testEnrichment(DM_l, lgS)
# Average 48 SNVs by kb of DM genes, average 18 SNVs by kb of non-DM genes. 

testEnrichment(DM_l, lgi)
# DM genes have on average 7 indels per kb, whereas DM have 4 indels per kb.

# In cortex
# load genes-variants overlap
cgS <- read.table("cortex_genesSNVsoverlap_BXDvsB6mm10permissive.txt", header=FALSE, stringsAsFactors=FALSE)
cgi <- read.table("cortex_genesindelsoverlap_BXDvsB6mm10permissive.txt", header=FALSE, stringsAsFactors=FALSE)

# load DM data
DM_c <- read.table("DifferentialMappingBXDvsB6mm10permissive_Limma_Cortex.txt", header=TRUE, stringsAsFactors=FALSE)

testEnrichment(DM_c, cgS)
# In cortex, average 43 SNVs by kb of DM genes, average 18 SNVs by kb of non-DM genes. 

testEnrichment(DM_c, cgi)
# In cortex, DM genes have on average 7 indels per kb, whereas non-DM genes have 4 indels per kb.

####################

# same for D2assembly vs B6mm10

# load genes-variants overlap
lgS <- read.table("liver_genesSNVsoverlap_D2assemblyvsB6mm10.txt", header=FALSE, stringsAsFactors=FALSE)
lgi <- read.table("liver_genesindelsoverlap_D2assemblyvsB6mm10.txt", header=FALSE, stringsAsFactors=FALSE)

# load DM data
DM_l <- read.table("DifferentialMapping_Limma_Liver_20190917.txt", header=TRUE, stringsAsFactors=FALSE)

testEnrichment(DM_l, lgS)
# Average 29 SNVs by kb of DM genes, average 15 SNVs by kb of non-DM genes. 

testEnrichment(DM_l, lgi)
# DM genes have on average 6 indels per kb, whereas DM have 3 indels per kb.

# In cortex
# load genes-variants overlap
cgS <- read.table("cortex_genesSNVsoverlap_D2assemblyvsB6mm10.txt", header=FALSE, stringsAsFactors=FALSE)
cgi <- read.table("cortex_genesindelsoverlap_D2assemblyvsB6mm10.txt", header=FALSE, stringsAsFactors=FALSE)

# load DM data
DM_c <- read.table("DifferentialMapping_Limma_Cortex_20190917.txt", header=TRUE, stringsAsFactors=FALSE)

testEnrichment(DM_c, cgS)
# In cortex, average 21 SNVs by kb of DM genes, average 16 SNVs by kb of non-DM genes. 

testEnrichment(DM_c, cgi)
# In cortex, DM genes have on average 4 indels per kb, whereas non-DM genes have 4 indels per kb.


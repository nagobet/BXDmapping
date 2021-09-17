# Retrieve statistics on D2-specific SNVs

# 30 March 2020

# prepare data tables
# Distance table
DistanceStats <- array(NA, c(21,2))
rownames(DistanceStats) <- c(1:19,"X","Y")
colnames(DistanceStats) <- c("NumberPairs","CloseVariant")
# Heterozygosity table
HeterozygosityStats <- array(NA, c(21,3))
rownames(HeterozygosityStats) <- c(1:19,"X","Y")
colnames(HeterozygosityStats) <- c("NumberVariants","HetVariant","HetClose")

# Define threshold for closeness between variants based on reads length.
closenessThreshold <- 100

# load data, calculate and store statistics
for(chr in c(1:19,"X","Y")){
  # load variants
  SNV <- read.table(paste0("F:/BXD/data/genome/D2specificVariants/", chr, ".tab"), stringsAsFactors=FALSE)
  # calculate distance
  SNVdistance <- SNV$V2[-length(SNV$V2)] - SNV$V2[-1]
  # store distance statistics
  DistanceStats[chr, "CloseVariant"] <- length(which(abs(SNVdistance)<closenessThreshold))
  DistanceStats[chr, "NumberPairs"] <- length(SNVdistance)
  # calculate heterozygosity
  Het <- grep("Het", SNV$V7)
  # count heterozygote SNV close to another heterozygote SNV
  SNVHetdistance <- SNV$V2[Het][-length(Het)] - SNV$V2[Het][-1]
  
  # store heterozygosity statistics
  HeterozygosityStats[chr, "HetVariant"] <- length(Het)
  HeterozygosityStats[chr, "NumberVariants"] <- length(SNV$V2)
  HeterozygosityStats[chr, "HetClose"] <- length(which(abs(SNVHetdistance)<closenessThreshold))
}

# save statisics as files (tab separeted format)
write.table(DistanceStats, file="F:/BXD/data/genome/D2specificVariants/DistanceStatsSNV.tsv", col.names=NA, sep="\t")
write.table(HeterozygosityStats, file="F:/BXD/data/genome/D2specificVariants/HeterozygosityStatsSNV.tsv", col.names=NA, sep="\t")

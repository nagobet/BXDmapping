# Extract D2 blocks from genotypes for BXD lines 

#GOAL: Extract D2 blocks from GN genotypes mm10 for all BXD lines 

# load function
source("F:/BXD/analysis/scripts/lowlayer/extractD2Blocks.R")

# BXD lines to 
##line <- "BXD43"
lines <- c("BXD43", "BXD44", "BXD45", "BXD48", "BXD49", 
           "BXD50", "BXD51", "BXD55", "BXD56", 
           "BXD61", "BXD64", "BXD65", "BXD66", "BXD67", 
           "BXD70", "BXD71", "BXD73", "BXD75", "BXD79", 
           "BXD81", "BXD83", "BXD84", "BXD85", "BXD87", "BXD89", 
           "BXD90", "BXD95", "BXD96", "BXD97", "BXD98", 
           "BXD100", "BXD101", "BXD103")

# extract and save in file
for(line in lines){
  write.table(extractD2Blocks(line), file=paste0("F:/BXD/data/genome/",line,"_D2blocks.bed"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}

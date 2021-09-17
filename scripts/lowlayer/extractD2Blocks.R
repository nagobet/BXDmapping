#GOAL: create R function to extract D2 blocks from genotypes in vcf

#DATE: 28 April 2020

# loading genotypes
genotypesmm10 <- read.table("F:/BXD/data/genome/GenotypesGNmm10.vcf", header=TRUE, comment.char="", stringsAsFactors=FALSE)
names(genotypesmm10)[names(genotypesmm10)=="BXD48a"] <- "BXD96"
names(genotypesmm10)[names(genotypesmm10)=="BXD65a"] <- "BXD97"
names(genotypesmm10)[names(genotypesmm10)=="BXD73b"] <- "BXD103"

##plot(genotypesmm10[genotypesmm10$X.CHROM==1,"POS"], genotypesmm10[genotypesmm10$X.CHROM==1,"X.CHROM"])
##plot(genotypesmm10[genotypesmm10$X.CHROM==1,"POS"], genotypesmm10[genotypesmm10$X.CHROM==1,"X.CHROM"], col=as.factor(genotypesmm10[genotypesmm10$X.CHROM==1,"BXD43"]))
##plot(genotypesmm10[genotypesmm10$X.CHROM==1,"POS"], genotypesmm10[genotypesmm10$X.CHROM==1,"X.CHROM"], col=as.factor(genotypesmm10[genotypesmm10$X.CHROM==1,"BXD98"]))

##blocks <- c("chrom", "chromStart", "chromEnd")
##line <- "BXD43"
##chr <- 1
##chr <- 2
## for(chr in levels(as.factor(genotypesmm10$X.CHROM))){
##   D2markers_idx <- grep("1/1",genotypesmm10[genotypesmm10$X.CHROM==chr,line])
##   ends <- D2markers_idx[c(D2markers_idx,NA)-(c(NA,D2markers_idx))!=1]
##   starts <- D2markers_idx[which(c(D2markers_idx,NA)-(c(NA,D2markers_idx))!=1)-1]
##   blocks_idx <- sort(c(range(D2markers_idx), starts, ends))
##   blocks_positions <- genotypesmm10$POS[genotypesmm10$X.CHROM==chr][blocks_idx]
##   blocks <- rbind(blocks, cbind(chr, matrix(blocks_positions, ncol=2, byrow=TRUE)))
## }
##blocks

extractD2Blocks <- function(line){
  # initiate object with header
  blocks <- c("chrom", "chromStart", "chromEnd")
  for(chr in c(1:19,"X")){
    # get index of start and end of D2 blocks
    D2markers_idx <- grep("1/1",genotypesmm10[genotypesmm10$X.CHROM==chr,line])
    ends <- D2markers_idx[c(D2markers_idx,NA)-(c(NA,D2markers_idx))!=1]
    starts <- D2markers_idx[which(c(D2markers_idx,NA)-(c(NA,D2markers_idx))!=1)-1]
    blocks_idx <- sort(c(range(D2markers_idx), starts, ends))
    # convert index into chromosomal positions
    blocks_positions <- genotypesmm10$POS[genotypesmm10$X.CHROM==chr][blocks_idx]
    blocks <- rbind(blocks, cbind(chr, matrix(blocks_positions, ncol=2, byrow=TRUE)))
    dimnames(blocks) <- list()
  }
  return(blocks)
}

# examples
##extractD2Blocks(line="BXD43")
##extractD2Blocks(line="BXD44")
##extractD2Blocks(line="BXD45")

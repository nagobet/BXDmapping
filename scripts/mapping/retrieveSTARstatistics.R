############################
# Retrieve STAR statistics #
############################

# 3 April 2020

retrieveStats <- function(ref, sample, path, statsName){
  # load
  statsdata <- read.table(paste0(path,"/", ref,"/",sample,"Log.final.out"), sep="\t", fill=TRUE, strip.white=TRUE, stringsAsFactors=FALSE)
  # format
  # remove unwanted rows
  data <- statsdata[-c(1:4,grep('^[A-Z -]+:$', statsdata$V1, value=FALSE, perl=TRUE, ignore.case=FALSE)),]
  # transform % into numeric
  data$V2 <- as.numeric(gsub("%", "",data$V2))
  # retrieve statistics of interest
  rownames(data) <- data$V1 
  return(data[statsName,2])
}

# Simple example
##path <- "F:/BXD/data/transcriptome/2_Mapping_STAR/OnPersonalizedRef/"
##statsName <- "Uniquely mapped reads % |"
##sample <- "DB1nsd"
##ref <- "DBA_2J_genome_nonrandomized_indels_maternal"
##retrieveStats(ref, paste0(sample, "_"), path, statsName)

# vectorize the function
retrieveStats <- Vectorize(FUN=retrieveStats, SIMPLIFY=TRUE, USE.NAMES=TRUE)
##retrieveStats(ref=dir(path, pattern="DBA_2J"), path="F:/BXD/data/transcriptome/2_Mapping_STAR/OnPersonalizedRef/", statsName="Uniquely mapped reads % |", sample=c("DB1nsd","LDB2sd"))

# vector example
##path <- "F:/BXD/data/transcriptome/2_Mapping_STAR/OnPersonalizedRef/"
##statsName <- "Uniquely mapped reads % |"
##references <- dir(path, pattern="DBA_2J_genome_")
##ref <- "DBA_2J_genome_nonrandomized_indels_maternal"
##samples <- gsub("_Log.final.out","", list.files(paste0(path, "/", ref), pattern="_Log.final.out"))

##res <- outer(FUN=retrieveStats, X=references, Y=paste0(samples, "_"), path="F:/BXD/data/transcriptome/2_Mapping_STAR/OnPersonalizedRef/", statsName="Uniquely mapped reads % |")
##res <- outer(FUN=retrieveStats, X=references, Y=paste(samples,"exactunique", sep="_"), path="F:/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/", statsName="% of reads unmapped: too many mismatches |")
##rownames(res) <- references
##colnames(res) <- samples
##res

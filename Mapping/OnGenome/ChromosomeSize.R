# size major chromosomes D2 vs B6 mm10

B6mm10_primaryassembly <- read.table("F:/BXD/references/genome/star2.7.0e_B6mm10_primaryassembly/chrNameLength.txt", stringsAsFactors=FALSE, row.names=1)
##B6mm10 <- read.table("F:/BXD/references/genome/star2.7.0e_B6mm10/chrNameLength.txt", stringsAsFactors=FALSE, row.names=1)
D2 <- read.table("F:/BXD/references/genome/star2.7.0e_D2/chrNameLength.txt", stringsAsFactors=FALSE, row.names=1)
sum(B6mm10_primaryassembly[as.character(1:19),"V2"])/10^9
##sum(B6mm10[as.character(1:19),"V2"])/10^9
sum(D2[as.character(1:19),"V2"])/10^9

B6mm10_primaryassembly[c(1:19, "X"),"V2"]-D2[c(1:19, "X"),"V2"]

barplot(t(cbind(B6mm10_primaryassembly[c(1:19, "X"),"V2"],D2[c(1:19, "X"),"V2"])), col=c("grey20","tan"), ylab="Size [bp]", xlab="Chromosome", names.arg=c(1:19, "X"), beside=TRUE)

# look at N bases

# count N bases faCount utility from UCSC https://www.biostars.org/p/19426/
#/home/ngobet/software/faCount references/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa > references/D2fastats.tab
#/home/ngobet/software/faCount references/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa > references/B6mm10primaryassemblyfastats.tab

B6mm10_primaryassembly <- read.table("F:/BXD/references/B6mm10primaryassemblyfastats.tab", stringsAsFactors=FALSE, header=TRUE, row.names=1, comment.char="")
D2 <- read.table("F:/BXD/references/D2fastats.tab", stringsAsFactors=FALSE, header=TRUE, row.names=1, comment.char="")

# full length
barplot(t(cbind(B6mm10_primaryassembly[c(1:19, "X"),"len"],D2[c(1:19, "X"),"len"])), col=c("grey20","tan"), ylab="Size [bp]", xlab="Chromosome", names.arg=c(1:19, "X"), beside=TRUE)
# non N
barplot(t(cbind(B6mm10_primaryassembly[c(1:19, "X"),"len"]-B6mm10_primaryassembly[c(1:19, "X"),"N"],D2[c(1:19, "X"),"len"]-D2[c(1:19, "X"),"N"])), col=c("grey20","tan"), ylab="Size [bp]", xlab="Chromosome", names.arg=c(1:19, "X"), beside=TRUE, main="size without N")
# N
barplot(t(cbind(B6mm10_primaryassembly[c(1:19, "X"),"N"],D2[c(1:19, "X"),"N"])), col=c("grey20","tan"), ylab="Size [bp]", xlab="Chromosome", names.arg=c(1:19, "X"), beside=TRUE, main="N")

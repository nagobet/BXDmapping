# PCA on mapping parameters

# DATE: 17 May 2020

# GOAL: understanding more important mapping parameters by doing a PCA on parameters for sample 67sd

# 1) loading
# load data parameters (explaining variables) and outcome (uniquely mapping reads percentage)
outcome <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/MappingParametersOptimization/valuesUniq.txt", stringsAsFactors=FALSE)
variables <- read.table("F:/BXD/data/transcriptome/2_mapping_STAR/MappingParametersOptimization/variables.txt", sep="_", stringsAsFactors=FALSE)
colnames(variables) <- c("annotation", "trimming", "deletion", "insertion", "intronMax", "sample", "mismatches", "Log")
outcome$V1 <- as.numeric(gsub("%", "", outcome$V1))
outcome_m <- cbind(outcome$V1[variables$sample=="L50nsd"], outcome$V1[variables$sample=="67sd"])
summary(outcome_m)
summary(variables)

# load packages
library(vioplot)
##install.packages('caroline')
library(caroline)

# 2) PCA analysis

PC <- prcomp(as.matrix(outcome_m))
summary(PC)

palette(c("grey", "red", "green3", "blue", "cyan", "magenta", "yellow"))
for(var in c("annotation", "trimming", "deletion", "insertion", "intronMax", "mismatches")){
  pairs(PC$x, pch=19, cex.lab=1.2, cex.main=2, cex.axis=1.2, col=as.factor(variables[,var]), main=var)
  plot(rownames(outcome), outcome$V1, pch=19, col=as.factor(variables[,var]), main=var)
  boxplot(outcome$V1 ~ as.factor(variables[,var]), ylab="Uniquely mapping reads %", col=1:length(levels(as.factor(variables[,"mismatches"]))), main=var, las=1)
  data <- list()
  for(l in levels(as.factor(variables[, var]))){
    data <- append(data, list(subset(outcome, variables[,var]==l)))
  }
  vioplot(data,col=1:length(levels(as.factor(variables[,var]))), main=var, las=1, ylab="Uniquely mapping reads %", names=levels(as.factor(variables[, var])))
}

# Importance of variables for outcome
fit1 <- aov(outcome$V1 ~ variables$annotation * variables$trimming * variables$deletion * variables$insertion * variables$intronMax * variables$mismatches)
summary(fit1)
TukeyHSD(fit1)
##plot(fit1)
fit2 <- aov(outcome$V1 ~ variables$annotation + variables$trimming + variables$deletion + variables$insertion + variables$intronMax + variables$mismatches)
summary(fit2)
TukeyHSD(fit2)

#CONCLUSION: It seems that annotation (with or without), trimming (Local or EndToEnd) and intronMax (0, 1) seems the more important variables.

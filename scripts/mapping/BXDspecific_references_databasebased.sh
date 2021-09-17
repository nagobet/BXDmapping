#!/bin/bash

#########################################################################
# BXD line-specific references based on D2 specific variants from dbSNP #
#########################################################################

# Get files from the Mouse Genome Project
# Store them in F:/BXD/data/genome/D2specificVariants

# Readme file for SNPs and indels
#ftp://ftp-mouse.sanger.ac.uk/current_snps/README
#renamed README_SNPs_Indels

# D2-specific SNPs, version 5
#ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz
# D2-specific indels, version 5
#ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf.gz

# all strains SVs
#ftp://ftp-mouse.sanger.ac.uk/current_svs/28strains.REL-1410-SV.sdp.tab.gz

# SVs info
#ftp://ftp-mouse.sanger.ac.uk/current_svs/README # renamed README_SVs
#ftp://ftp-mouse.sanger.ac.uk/current_svs/liftOver_definitions
#ftp://ftp-mouse.sanger.ac.uk/current_svs/svs.mm9.unmapped_mm10.txt

#ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd118_Keane_et_al_2011/gvf/estd118_Keane_et_al_2011.2013-01-07.MGSCv37.gvf

# 11 March 2020

# GOAL: test different personalized references with D2-specific SNPs and indels from dbSNP

# define frequently used directories
refdir=/mnt/nas/BXD/data/PersonalizedReferences
aligndir=/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnPersonalizedRef
tmpdir=/mnt/md0/BXD/TMP

line=DBA_2J
# echo $line
# # build personalized genomes
# export JAVA_HOME=/home/ngobet/software/java/jdk-13.0.2
# export PATH=$JAVA_HOME/bin:$PATH
# # unphased variant randomized or not
# mkdir $refdir/$line\_randomized_indels $refdir/$line\_nonrandomized_indels
# cd $refdir/$line\_randomized_indels
# java -jar /home/ngobet/software/vcf2diploid-master/vcf2diploid.jar \
# -id $line \
# -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
# -vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf  > vcf2diploid.out 2> vcf2diploid.err
# cd $refdir/$line\_nonrandomized_indels
# java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
# -id $line \
# -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
# -vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf > vcf2diploid.out 2> vcf2diploid.err
# mkdir $refdir/$line\_randomized_indelsSNVs $refdir/$line\_nonrandomized_indelsSNVs
# cd $refdir/$line\_randomized_indelsSNVs
# java -jar /home/ngobet/software/vcf2diploid-master/vcf2diploid.jar \
# -id $line \
# -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
# -vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > vcf2diploid.out 2> vcf2diploid.err
# cd $refdir/$line\_nonrandomized_indelsSNVs
# java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
# -id $line \
# -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
# -vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > vcf2diploid.out 2> vcf2diploid.err
# mkdir $refdir/$line\_randomized_SNVs $refdir/$line\_nonrandomized_SNVs
# cd $refdir/$line\_randomized_SNVs
# java -jar /home/ngobet/software/vcf2diploid-master/vcf2diploid.jar \
# -id $line \
# -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
# -vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > vcf2diploid.out 2> vcf2diploid.err
# cd $refdir/$line\_nonrandomized_SNVs
# java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
# -id $line \
# -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
# -vcf /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > vcf2diploid.out 2> vcf2diploid.err

# download unplaced chromosomes
#ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa.gz

# build references index for STAR
for parental in maternal paternal
do
	for unphased in nonrandomized randomized
	do
	for variants in SNVs indelsSNVs indels
	do
	echo $parental\_$unphased
	mkdir $refdir/$line\_$unphased\_$variants/star_$parental $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation
	# liftover transcriptome annotation
	##/home/ngobet/software/liftOver -gff /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf $refdir/$line\_$unphased\_$variants/$parental\.chain $refdir/$line\_$unphased\_$variants/$parental\.gtf $refdir/$line\_$unphased\_$variants/$parental\_unlifted.gtf
	# make genome reference STAR index
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental/star_$parental\_ > $refdir/$line\_$unphased\_$variants/star_$parental/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental/indexing.err
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --sjdbGTFfile $refdir/$line\_$unphased\_$variants/$parental.gtf --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/star_$parental\_withannotation_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.err
	# Run STAR alignment
	mkdir $tmpdir/$line\_genome_$unphased\_$variants\_$parental/ $tmpdir/$line\_genomeandann_$unphased\_$variants\_$parental/
	##for sample in DB1nsd DB1sd DB2nsd DB2sd LDB1nsd LDB1sd LDB2nsd LDB2sd
	for sample in DB1nsd DB2nsd DBnsd BDnsd B61nsd B62nsd LDB1nsd LDB2nsd LDBnsd LBDnsd LB61nsd LB62nsd
	do
	##sample=L98nsd
	# restrictive (exact matches) on genome
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental --outFileNamePrefix $tmpdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 > $tmpdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_align_restrictive.out 2> $tmpdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_align_restrictive.err
	# defaults on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation --outFileNamePrefix $tmpdir/$line\_genomeandann_$unphased\_$variants\_$parental/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --quantMode GeneCounts > $tmpdir/$line\_genomeandann_$unphased\_$variants\_$parental/$sample\_align_default.out 2> $tmpdir/$line\_genomeandann_$unphased\_$variants\_$parental/$sample\_align_default.err
	done
	# move file from temporary to output directory
	##mv $tmpdir/$line\_genome_$unphased\_$variants\_$parental $aligndir/
	mv $tmpdir/$line\_genomeandann_$unphased\_$variants\_$parental/ $aligndir/
	done
	done
done

# STAR alignment on maternal or paternal, randomized or not.
# STAR alignment parameters restrictive 
# different samples: B6, D2, 2 F1, a few BXD: SD and NSD, cortex and liver

# use kallisto?

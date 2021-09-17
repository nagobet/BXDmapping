#!/bin/bash

####################################################
# BXD line-specific references based on genotyping #
####################################################

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

# GOAL: test different personalized references using GeneNetwork genotypes

# define frequently used directories
refdir=/mnt/nas/BXD/data/PersonalizedReferences
aligndir=/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnPersonalizedRef
tmpdir=/mnt/md0/BXD/TMP

##for line in DBA_2J BXD43 BXD48 BXD49 BXD83 BXD98 BXD100 B6D2F1 D2B5F1
##for line in DBA_2J B6D2F1 BXD43 BXD48 BXD49 BXD83 BXD98 BXD100 
# strains chosen from phylogeny https://ars.els-cdn.com/content/image/1-s2.0-S2405471217304866-gr1_lrg.jpg
#all parental + F1
#BXD48 BXD98 are closer to B6
#BXD43 BXD100 are closer to D2
#BXD49 BXD83 are intermediate
##lines=($(cut -f 6 data/ConvertLineNames.tsv | grep "BXD" | awk 'NR > 4 { print }'))
##for line in "${lines[@]}"
for line in B6D2F1
do
	echo $line
	##linenumber=$(grep -o [0-9]* <<< $line)
	##echo $linenumber
	# build personalized genomes
	##export JAVA_HOME=/home/ngobet/software/java/jdk-13.0.2
	##export PATH=$JAVA_HOME/bin:$PATH
	# unphased variant randomized or not
	##mkdir $refdir/$line\_randomized_genotypes $refdir/$line\_nonrandomized_genotypes
	##cd $refdir/$line\_randomized_genotypes
	##java -jar /home/ngobet/software/vcf2diploid-master/vcf2diploid.jar \
	## -id $line \
	## -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	## -vcf /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf > vcf2diploid.out 2> vcf2diploid.err
	## cd $refdir/$line\_nonrandomized_genotypes
	## java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
	## -id $line \
	## -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	## -vcf /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf > vcf2diploid.out 2> vcf2diploid.err
	##for parental in maternal paternal
	for parental in paternal
	do
	##for unphased in nonrandomized randomized
	for unphased in nonrandomized
	do
	echo $parental\_$unphased\_genotypes
	##mkdir $refdir/$line\_$unphased\_genotypes/star_$parental $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation
	# liftover transcriptome annotation
	##/home/ngobet/software/liftOver -gff /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf $refdir/$line\_$unphased\_genotypes/$parental\.chain $refdir/$line\_$unphased\_genotypes/$parental\.gtf $refdir/$line\_$unphased\_genotypes/$parental\_unlifted.gtf
	# build references index for STAR
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_genotypes/star_$parental --genomeFastaFiles $refdir/$line\_$unphased\_genotypes/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --outFileNamePrefix $refdir/$line\_$unphased\_genotypes/star_$parental/star_$parental\_ > $refdir/$line\_$unphased\_genotypes/star_$parental/indexing.out 2> $refdir/$line\_$unphased\_genotypes/star_$parental/indexing.err
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation --genomeFastaFiles $refdir/$line\_$unphased\_genotypes/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --sjdbGTFfile $refdir/$line\_$unphased\_genotypes/$parental.gtf --outFileNamePrefix $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation/star_$parental\_withannotation_ > $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation/indexing.out 2> $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation/indexing.err
	# Run STAR alignment
	mkdir $tmpdir/$line\_genome_$unphased\_genotypes\_$parental/ $tmpdir/$line\_genomeandann_$unphased\_genotypes\_$parental/
	##for sample in 43nsd 43sd 48nsd 48sd 49nsd 49sd 83nsd 83sd 98nsd 98sd 100nsd 100sd B61nsd B61sd B62nsd B62sd BDnsd BDsd DB1nsd DB1sd DB2nsd DB2sd DBnsd DBsd L43nsd L43sd L48nsd L48sd L49nsd L49sd L83nsd L83sd L98nsd L98sd L100nsd L100sd LB61nsd LB61sd LB62nsd LB62sd LBDnsd LBDsd LDB1nsd LDB1sd LDB2nsd LDB2sd LDBnsd LDBsd
	##for sample in DB1nsd DB2nsd DBnsd BDnsd B61nsd B62nsd LDB1nsd LDB2nsd LDBnsd LBDnsd LB61nsd LB62nsd
	##samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*$linenumber*.fastq.gz | grep -o "[L]\{0,1\}[0-9]\{2,3\}[nsd]\{2,3\}*"))
	samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L]\{0,1\}[0-9BDnsd]\{3,6\}*"))
	for sample in "${samples[@]}"
	do
	echo $sample
	##sample=L98nsd
	# restrictive (exact matches) on genome
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_genotypes/star_$parental --outFileNamePrefix $tmpdir/$line\_genome_$unphased\_genotypes\_$parental/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 > $tmpdir/$line\_genome_$unphased\_genotypes\_$parental/$sample\_align_restrictive.out 2> $tmpdir/$line\_genome_$unphased\_genotypes\_$parental/$sample\_align_restrictive.err
	# defaults on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation --outFileNamePrefix $tmpdir/$line\_genomeandann_$unphased\_genotypes\_$parental/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outStd Log --quantMode GeneCounts > $tmpdir/$line\_genomeandann_$unphased\_genotypes\_$parental/$sample\_align_default.out 2> $tmpdir/$line\_genomeandann_$unphased\_genotypes\_$parental/$sample\_align_default.err
	done
	# move file from temporary to output directory
	##mv $tmpdir/$line\_genome_$unphased\_genotypes\_$parental $aligndir/
	mv $tmpdir/$line\_genomeandann_$unphased\_genotypes\_$parental/ $aligndir/
	done
	done
	cd /mnt/nas/BXD/
done

# STAR alignment on maternal or paternal, randomized or not.
# STAR alignment parameters restrictive 
# different samples: B6, D2, 2 F1, a few BXD: SD and NSD, cortex and liver

# use kallisto?

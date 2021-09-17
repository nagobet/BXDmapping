#!/bin/bash

# Genotypes Imputation D2blockmethod

# GOAL: impute genotypes with my D2 blocks homemade simple method

# DATE: 19 November 2020

# Resource files
#data/genome/GenotypesGNmm10.vcf
#data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf
#data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf

# data/genome/BXD43_D2blocks.bed

# Mapping
refdir=/mnt/nas/BXD/data/PersonalizedReferences
outdir=/mnt/md0/BXD

##line=BXD49
parental=paternal
unphased=nonrandomized
variants=genotypesandimputed

# all lines
cut -f 6 data/ConvertLineNames.tsv | grep "BXD" | awk 'NR > 4 { print }' > data/lines.txt
for line in $(cat data/lines.txt)
do
	echo $line

	# Impute genotypes for BXD line
	# retrieve header
	grep "#CHROM" data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# adapt header to BXD line
	sed -i "s/DBA_2J/$line/g" data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# impute genotypes for BXD line
	grep -v "NA" data/genome/$line\_D2blocks.bed > blocks.bed
	/home/ngobet/software/bedtools2/bin/intersectBed -a data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf -b blocks.bed -wa >> data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	/home/ngobet/software/bedtools2/bin/intersectBed -a data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf -b blocks.bed -wa >> data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	rm blocks.bed
	
	# Build personalized genomes
	export JAVA_HOME=/home/ngobet/software/java/jdk-13.0.2
	export PATH=$JAVA_HOME/bin:$PATH
	# unphased variant randomized or not
	mkdir $refdir/$line\_$unphased\_$variants
	cd $refdir/$line\_$unphased\_$variants
	java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
	 -id $line \
	 -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	 -vcf /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf /mnt/nas/BXD/data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf > $refdir/$line\_$unphased\_$variants/vcf2diploid.out 2> $refdir/$line\_$unphased\_$variants/vcf2diploid.err
	cd /mnt/nas/BXD/
	mkdir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation
	# liftover transcriptome annotation
	/home/ngobet/software/liftOver -gff /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf $refdir/$line\_$unphased\_$variants/$parental\.chain $refdir/$line\_$unphased\_$variants/$parental\.gtf $refdir/$line\_$unphased\_$variants/$parental\_unlifted.gtf
	# build references index for STAR
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 16 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/star_$parental\_withoutannotation_\_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/indexing.err
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 4 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --sjdbGTFfile $refdir/$line\_$unphased\_$variants/$parental.gtf --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/star_$parental\_withannotation_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.err
	
	# # Run STAR alignment
	# mkdir $outdir/$line\_genome_$unphased\_$variants\_$parental/ $outdir/$line\_genomeandann_$unphased\_$variants\_$parental/
	# ##for sample in 43nsd 43sd L43nsd L43sd
	# linenumber=$(grep -o [0-9]* <<< $line)
	# samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDtnsd]\{4,7\}*" | grep "$linenumber" | grep -v "B"))
	# for sample in "${samples[@]}"
	# do
		# echo $sample
		# # restrictive (exact matches) on genome
		# /software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation --outFileNamePrefix $outdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 > $outdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_align_restrictive.out 2> $outdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_align_restrictive.err
		# # defaults on genome + annotation
		# /software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation --outFileNamePrefix $outdir/$line\_genomeandann_$unphased\_$variants\_$parental/$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outStd Log --quantMode GeneCounts > $outdir/$line\_genomeandann_$unphased\_$variants\_$parental/$sample\_align_default.out 2> $outdir/$line\_genomeandann_$unphased\_$variants\_$parental/$sample\_align_default.err
		# # remove unnecessary files
		# rm $outdir/$line\_genome*_$unphased\_$variants\_$parental/$sample\_*.bam
		# rm $outdir/$line\_genome*_$unphased\_$variants\_$parental/$sample\_SJ.out.tab
	# done
done

#!/bin/bash

# STAR mapping optimization

# GOAL: run mapping various STAR mapping parameters to identify optimal parameters

# DATE: 17 May 2020 - 26 June 2020

# INFO:
# Steps from reads of samples (.fastq.gz) to cis eQTL detection
# Some steps are run on Wally cluster (see MappingOptimization_mainOnCluster.sh)

###################################################################################

# Part 1: Mapping

# Constant definitions
refdir=/home/ngobet/projects/BXD/references
outdir=/mnt/md0/BXD/MappingOptimization
parental=paternal
unphased=nonrandomized
variants=genotypesandimputed
deletion=2
insertion=2

# prepare mapping commands
mkdir $outdir
for line in $(cat data/lines.txt)
do
	echo "#$line"
	echo "#$line" > /mnt/nas/BXD/analysis/scripts/mapping/MappingOptimization_$line\_listcommands.sh
	echo "#DATE:" `date` >> /mnt/nas/BXD/analysis/scripts/mapping/MappingOptimization_$line\_listcommands.sh
	# Run STAR alignment
	for annotation in withoutannotation withannotation
	do
		echo "#$annotation $count"
		for trimming in Local EndToEnd
		do
			echo "##trimming: $trimming"
			for intronMax in 1 0
			do
				echo "###intronMax: $intronMax"
				for mismatches in 0 1 2 3 10
				do
					echo "####mismatches: $mismatches"
					mkdir -p $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/
					linenumber=$(grep -o [0-9]* <<< $line)
					samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9D]\{2,4\}[n]\{0,1\}sd*" | grep "$linenumber" | sort | uniq))
					for sample in "${samples[@]}"
					do
						# mapping command (alignment with STAR + gene counting)
						# mapping without transcriptome annotation
						if [ $annotation = "withoutannotation" ]
						then
							echo "/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_$annotation --outFileNamePrefix $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_ --outSAMmode Full --outStd SAM --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --scoreDelOpen -$deletion --scoreInsOpen -$insertion --alignIntronMax $intronMax --alignEndsType $trimming --outFilterMismatchNmax $mismatches 2> $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_align.err | htseq-count -f bam -q -s reverse -t exon -m union - $refdir/$line\_$unphased\_$variants/$parental.gtf 2> $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_htseqcount.err 1> $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_ReadsPerGene.out.tab; rm $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_SJ.out.tab" >> /mnt/nas/BXD/analysis/scripts/mapping/MappingOptimization_$line\_listcommands.sh
						# mapping with transcriptome annotation
						elif [ $annotation = "withannotation" ]
						then
							echo "/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_$annotation --outFileNamePrefix $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_ --outSAMtype BAM Unsorted --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --quantMode GeneCounts --scoreDelOpen -$deletion --scoreInsOpen -$insertion --alignIntronMax $intronMax --alignEndsType $trimming --outFilterMismatchNmax $mismatches > $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_align_default.out 2> $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_align.err; rm $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_*.bam; rm $outdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_SJ.out.tab" >> /mnt/nas/BXD/analysis/scripts/mapping/MappingOptimization_$line\_listcommands.sh
						fi
					done #end samples loop
				done #end mismatches loop
			done #end intronMax loop
		done #end trimming loop
	done #end annotation loop
done #end lines loop

# Run mapping commands
for line in $(cat data/lines.txt)
##for line in BXD67 BXD50
do
	echo "#$line"
	# # Impute genotypes for BXD line
	# # retrieve header
	# grep "#CHROM" data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# # adapt header to BXD line
	# sed -i "s/DBA_2J/$line/g" data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# # impute genotypes for BXD line
	# /home/ngobet/software/bedtools2/bin/intersectBed -a data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf -b data/genome/$line\_D2blocks.bed -wa >> data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# /home/ngobet/software/bedtools2/bin/intersectBed -a data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf -b data/genome/$line\_D2blocks.bed -wa >> data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf

	# # Build personalized genomes
	# export JAVA_HOME=/home/ngobet/software/java/jdk-13.0.2
	# export PATH=$JAVA_HOME/bin:$PATH
	# # unphased variant randomized or not
	# mkdir $refdir/$line\_$unphased\_$variants
	# cd $refdir/$line\_$unphased\_$variants
	# java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
	 # -id $line \
	 # -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	 # -vcf /mnt/nas/BXD/data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf > $refdir/$line\_$unphased\_$variants/vcf2diploid.out 2> $refdir/$line\_$unphased\_$variants/vcf2diploid.err
	# cd /mnt/nas/BXD/
	# mkdir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation
	# # liftover transcriptome annotation
	# /home/ngobet/software/liftOver -gff /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf $refdir/$line\_$unphased\_$variants/$parental\.chain $refdir/$line\_$unphased\_$variants/$parental\.gtf $refdir/$line\_$unphased\_$variants/$parental\_unlifted.gtf
	# # build references index for STAR
	# /software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 16 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/star_$parental\_withoutannotation_\_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/indexing.err
	# /software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 16 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --sjdbGTFfile $refdir/$line\_$unphased\_$variants/$parental.gtf --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/star_$parental\_withannotation_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.err
	# Get gene lengths
	# GTFtools only consider genes on conventional genes (help from: https://genomespot.blogspot.com/2019/01/using-gtf-tools-to-get-gene-lengths.html)
	grep -v '#' data/PersonalizedReferences/$line\_nonrandomized_genotypesandimputed/paternal.gtf | awk '{OFS="\t"} $1=1' > tmp.gtf
	/home/ngobet/software/GTFtools_0.6.5/gtftools.py -l data/PersonalizedReferences/$line\_nonrandomized_genotypesandimputed/genelength.txt tmp.gtf
	rm tmp.gtf
	
	# copy reference directory to SSD memory
	cp -r data/PersonalizedReferences/$line\_$unphased\_$variants $refdir
	
	# run jobs
	nohup parallel --delay 1 -j 4 < /mnt/nas/BXD/analysis/scripts/mapping/MappingOptimization_$line\_listcommands.sh 1> $outdir/MappingOptimization_$line\_listcommands.out 2> $outdir/MappingOptimization_$line\_listcommands.err
	
	# remove reference directory copy on SSD memory
	rm -r $refdir/$line\_$unphased\_$variants
done

# # retrieve stats
# mkdir data/transcriptome/2_mapping_STAR/MappingParametersOptimization
# grep "reads %" /mnt/md0/BXD/BXD[56][07]_genotypesandimputed_*/*_Log.final.out | cut -f 2 > data/transcriptome/2_mapping_STAR/MappingParametersOptimization/valuesUniq.txt
# grep "reads %" /mnt/md0/BXD/BXD[56][07]_genotypesandimputed_*/*_Log.final.out | cut -d "_" -f 3-9 | cut -d "." -f 1 | sed 's/\//\_/g' > data/transcriptome/2_mapping_STAR/MappingParametersOptimization/variables.txt

module load R/3.4.2
listsettings=($(ls -d data/transcriptome/2_mapping_STAR/MappingParametersOptimization/*/))
for setting in "${listsettings[@]}"
do
	##echo setting is: $setting
	# group counts
	analysis/scripts/mapping/MergeCount.sh $setting _ReadsPerGene.out.tab

	# normalize counts CPM
	Rscript analysis/scripts/mapping/normalizeGeneCountsCPM.R $setting
	# normalize counts TPM
	Rscript analysis/scripts/mapping/normalizeGeneCountsTPM.R $setting
done

# prepare gene annotation for eQTL analysis
Rscript analysis/scripts/mapping/retrieveGenePosition.R

# Part 2: eQTL analysis [on cluster]

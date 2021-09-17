# Align RNA on transcriptome
# Using STAR
# exactly mapped, 1 mismatch allowed, 2, 3, 4, 5.

# 1 July 2019

# Real alignments

# set working directory
cd /mnt/nas/BXD

# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP/OnTranscriptome

# 2.1 Generating genome indexes files
mkdir -p references/star2.7.0e_transcriptome_B6_mm10
mkdir -p references/star2.7.0e_transcriptome_D2
# [waiting]

# generate transcriptome indexes for star 2.7.0.e
module add Alignment/STAR/2.7.0e; STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/star2.7.0e_transcriptome_B6_mm10/ --genomeFastaFiles references/Mus_musculus.GRCm38.dna.toplevel.fa --sjdbGTFfile references/Mus_musculus.GRCm38.94.gtf --outFileNamePrefix references/star2.7.0e_transcriptome_B6_mm10/B6_mm10_ &
##module add Alignment/STAR/2.7.0e; STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/star2.7.0e_transcriptome_D2/ --genomeFastaFiles references/Mus_musculus_dba2j.DBA_2J_v1.dna.toplevel.fa --sjdbGTFfile references/Mus_musculus_dba2j.DBA_2J_v1.94.gtf --outFileNamePrefix references/star2.7.0e_transcriptome_D2/D2_ &
# test alignement on transcriptome with strict parameters
module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/B61nsd\_exactunique --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/B61nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/align\_B61nsd\_exactunique.out 2> /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/align\_B61nsd\_exactunique.err &
module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/B61nsd\_exactuniqueIntrons --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/B61nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/align\_B61nsd\_exactuniqueIntrons.out 2> /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/align\_B61nsd\_exactuniqueIntrons.err & 
### not working, there are junctions (so introns) in the 2 cases. --sjdbOverhang parameter seems redefined by the reference genome used. => putting it back to 0:
module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/B61nsd\_exactunique2 --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/B61nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 --sjdbOverhang 0 > /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/align\_B61nsd\_exactunique2.out 2> /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/align\_B61nsd\_exactunique2.err &
### not working!!! (different --sjdbOverhang value for reference index generation and reads alignment is not accepted)



# generate list of onliner commands
for ref in D2
do
	# for all samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz; do
	# take only the parental and F1 samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*.fastq.gz; do
	# for BXD samples (every sample except parental and F1)
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/[^BD][^BD]*.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $fqfile
	##echo $sample
	spe="exactunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##spe="1mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##spe="2mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##spe="3mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##spe="4mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##spe="5mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/star2.7.0e_transcriptome_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommandsB6D2F1_alignRNAontranscriptome\_$spe.sh
	done
done

# 19 December 2019

# GOAL: align with STAR on the transcriptome sequences (to compare to Kallisto)

# Generate STAR indexes for transcriptome sequences
# create directories
mkdir /mnt/nas/BXD/references/transcriptome/STAR_2.7.0e_B6mm10
mkdir /mnt/nas/BXD/references/transcriptome/STAR_2.7.0e_D2
# uncompress fasta files
gzip -dk /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.cdna.all.fa.gz
gzip -dk /mnt/nas/BXD/references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.cdna.all.fa.gz
# generate index B6 mm10
nohup /software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 4 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir /mnt/nas/BXD/references/transcriptome/STAR_2.7.0e_B6mm10/ --genomeFastaFiles /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.cdna.all.fa --outFileNamePrefix references/transcriptome/STAR_2.7.0e_B6mm10/B6mm10_ &
#I had to increase the --limitGenomeGenerateRAM otherwise it crashed.
# generate index D2
nohup /software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 4 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir /mnt/nas/BXD/references/transcriptome/STAR_2.7.0e_D2/ --genomeFastaFiles /mnt/nas/BXD/references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.cdna.all.fa --outFileNamePrefix references/transcriptome/STAR_2.7.0e_D2/D2_ > references/transcriptome/STAR_2.7.0e_D2/indexingD2.out 2> references/transcriptome/STAR_2.7.0e_D2/indexingD2.err&

# Align reads to transcriptome only, allowing different number of mismatches

# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP
# generate list of onliner commands
for ref in D2 B6mm10
do
	mkdir $outdir/$ref
	mkdir /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref
	# for all samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz; do
	# take only the parental and F1 samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*.fastq.gz; do
	# for BXD samples (every sample except parental and F1)
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/[^BD][^BD]*.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $fqfile
	##echo $sample
	for m in {0..10}; do
	spe=$m"mismatch"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/transcriptome/STAR_2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax $m > $outdir/$ref/$sample\_$spe\_align.out 2> $outdir/$ref/$sample\_$spe\_align.err " >> scripts/listcommands_alignRNAontranscriptome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/listcommands_alignRNAontranscriptome\_$spe.sh
	done # end of loop for number of mismatches
	done # end of loop for samples
done # end of loop for reference transcriptome

# running parental and F1 alignments
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_0mismatch.sh 1> $outdir/0mismatch_parallel.out 2> $outdir/0mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_1mismatch.sh 1> $outdir/1mismatch_parallel.out 2> $outdir/1mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_2mismatch.sh 1> $outdir/2mismatch_parallel.out 2> $outdir/2mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_3mismatch.sh 1> $outdir/3mismatch_parallel.out 2> $outdir/3mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_4mismatch.sh 1> $outdir/4mismatch_parallel.out 2> $outdir/4mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_5mismatch.sh 1> $outdir/5mismatch_parallel.out 2> $outdir/5mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_10mismatch.sh 1> $outdir/10mismatch_parallel.out 2> $outdir/10mismatch_parallel.err &

nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_6mismatch.sh 1> $outdir/6mismatch_parallel.out 2> $outdir/6mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_7mismatch.sh 1> $outdir/7mismatch_parallel.out 2> $outdir/7mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_8mismatch.sh 1> $outdir/8mismatch_parallel.out 2> $outdir/8mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAontranscriptome_9mismatch.sh 1> $outdir/9mismatch_parallel.out 2> $outdir/9mismatch_parallel.err &


# 13 January 2020

# Align with STAR on transcriptome with default options

# generate commands
# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP
# generate list of onliner commands
for ref in D2 B6mm10
do
	mkdir $outdir/$ref
	mkdir /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref
	# for all samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	spe="default"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/transcriptome/STAR_2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 > $outdir/$ref/$sample\_$spe\_align.out 2> $outdir/$ref/$sample\_$spe\_align.err " >> scripts/mapping/listcommands_alignRNAontranscriptome\_$spe.sh
	echo "rm /home/ngobet/projects/BXD/data/TMP/$ref/*$spe*.bam" >> scripts/mapping/listcommands_alignRNAontranscriptome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnTranscriptome/$ref" >> scripts/mapping/listcommands_alignRNAontranscriptome\_$spe.sh
	done # end of loop for samples
done # end of loop for reference transcriptome

# 23 January 2020

# GOAL: retrieve mapping stats (uniquely and multi mapped reads) aligned with STAR on the transcriptome, for 0 to 10 mismatches (restrictive conditions) or default (permissive).

# Retrieve info for all samples
# prepare headers
echo "Sample_Name" > tmp_headers
echo "Total" >> tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
echo "Uniq6" >> tmp_headers
echo "Uniq7" >> tmp_headers
echo "Uniq8" >> tmp_headers
echo "Uniq9" >> tmp_headers
echo "Uniq10" >> tmp_headers
echo "Multi0" >> tmp_headers
echo "Multi1" >> tmp_headers
echo "Multi2" >> tmp_headers
echo "Multi3" >> tmp_headers
echo "Multi4" >> tmp_headers
echo "Multi5" >> tmp_headers
echo "Multi6" >> tmp_headers
echo "Multi7" >> tmp_headers
echo "Multi8" >> tmp_headers
echo "Multi9" >> tmp_headers
echo "Multi10" >> tmp_headers
echo "UniqDefault" >> tmp_headers
echo "MultiDefault" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on D2
grep "Number of input reads" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_0mismatchLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Number of input reads" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_0mismatchLog.final.out | cut -f 2 > tmp_total
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_0mismatchLog.final.out | cut -f 2 > tmp_uniq0
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_1mismatchLog.final.out | cut -f 2 > tmp_uniq1
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_2mismatchLog.final.out | cut -f 2 > tmp_uniq2
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_3mismatchLog.final.out | cut -f 2 > tmp_uniq3
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_4mismatchLog.final.out | cut -f 2 > tmp_uniq4
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_5mismatchLog.final.out | cut -f 2 > tmp_uniq5
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_6mismatchLog.final.out | cut -f 2 > tmp_uniq6
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_7mismatchLog.final.out | cut -f 2 > tmp_uniq7
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_8mismatchLog.final.out | cut -f 2 > tmp_uniq8
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_9mismatchLog.final.out | cut -f 2 > tmp_uniq9
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_10mismatchLog.final.out | cut -f 2 > tmp_uniq10
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_defaultLog.final.out | cut -f 2 > tmp_uniqdefault
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_0mismatchLog.final.out | cut -f 2 > tmp_multi0
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_1mismatchLog.final.out | cut -f 2 > tmp_multi1
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_2mismatchLog.final.out | cut -f 2 > tmp_multi2
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_3mismatchLog.final.out | cut -f 2 > tmp_multi3
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_4mismatchLog.final.out | cut -f 2 > tmp_multi4
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_5mismatchLog.final.out | cut -f 2 > tmp_multi5
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_6mismatchLog.final.out | cut -f 2 > tmp_multi6
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_7mismatchLog.final.out | cut -f 2 > tmp_multi7
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_8mismatchLog.final.out | cut -f 2 > tmp_multi8
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_9mismatchLog.final.out | cut -f 2 > tmp_multi9
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_10mismatchLog.final.out | cut -f 2 > tmp_multi10
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/D2/*_defaultLog.final.out | cut -f 2 > tmp_multidefault
paste tmp_sample tmp_total tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq6 tmp_uniq7 tmp_uniq8 tmp_uniq9 tmp_uniq10 tmp_multi0 tmp_multi1 tmp_multi2 tmp_multi3 tmp_multi4 tmp_multi5 tmp_multi6 tmp_multi7 tmp_multi8 tmp_multi9 tmp_multi10 tmp_uniqdefault tmp_multidefault > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnTranscriptome/MappingStatisticsD2.tsv
rm tmp_*

# retrieve statistics for alignment on B6 mm10
grep "Number of input reads" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_0mismatchLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Number of input reads" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_0mismatchLog.final.out | cut -f 2 > tmp_total
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_0mismatchLog.final.out | cut -f 2 > tmp_uniq0
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_1mismatchLog.final.out | cut -f 2 > tmp_uniq1
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_2mismatchLog.final.out | cut -f 2 > tmp_uniq2
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_3mismatchLog.final.out | cut -f 2 > tmp_uniq3
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_4mismatchLog.final.out | cut -f 2 > tmp_uniq4
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_5mismatchLog.final.out | cut -f 2 > tmp_uniq5
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_6mismatchLog.final.out | cut -f 2 > tmp_uniq6
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_7mismatchLog.final.out | cut -f 2 > tmp_uniq7
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_8mismatchLog.final.out | cut -f 2 > tmp_uniq8
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_9mismatchLog.final.out | cut -f 2 > tmp_uniq9
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_10mismatchLog.final.out | cut -f 2 > tmp_uniq10
grep "Uniquely mapped reads number" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_defaultLog.final.out | cut -f 2 > tmp_uniqdefault
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_0mismatchLog.final.out | cut -f 2 > tmp_multi0
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_1mismatchLog.final.out | cut -f 2 > tmp_multi1
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_2mismatchLog.final.out | cut -f 2 > tmp_multi2
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_3mismatchLog.final.out | cut -f 2 > tmp_multi3
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_4mismatchLog.final.out | cut -f 2 > tmp_multi4
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_5mismatchLog.final.out | cut -f 2 > tmp_multi5
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_6mismatchLog.final.out | cut -f 2 > tmp_multi6
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_7mismatchLog.final.out | cut -f 2 > tmp_multi7
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_8mismatchLog.final.out | cut -f 2 > tmp_multi8
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_9mismatchLog.final.out | cut -f 2 > tmp_multi9
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_10mismatchLog.final.out | cut -f 2 > tmp_multi10
grep "Number of reads mapped to multiple loci" data/transcriptome/2_Mapping_STAR/OnTranscriptome/B6mm10/*_defaultLog.final.out | cut -f 2 > tmp_multidefault
paste tmp_sample tmp_total tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq6 tmp_uniq7 tmp_uniq8 tmp_uniq9 tmp_uniq10 tmp_multi0 tmp_multi1 tmp_multi2 tmp_multi3 tmp_multi4 tmp_multi5 tmp_multi6 tmp_multi7 tmp_multi8 tmp_multi9 tmp_multi10 tmp_uniqdefault tmp_multidefault > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnTranscriptome/MappingStatisticsB6mm10.tsv
rm tmp_*

rm headers

# info on STAR input

# unmapped others = STAR cannot find anchors (= exactly matching seed) to map at not too many places. Likely contamination with other species, low complexity loci, or very poor sequencing quality.
# source: https://groups.google.com/forum/#!topic/rna-star/Bo_m61aOabs
# % of reads unmapped: too short = alignments too shorts. Could be caused by over-trimmed reads (read length too short) or normal read not mapping well.
# source: https://github.com/alexdobin/STAR/issues/164
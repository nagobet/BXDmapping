# align RNA on genome

# 12 June 2019

# set working directory
cd /mnt/nas/BXD/

# Generate the reference genomes (without annotation) for star

# create directories
mkdir -p references/genome/star2.7.0e_B6mm10
mkdir -p references/genome/star2.7.0e_B6mm9
mkdir -p references/genome/star2.7.0e_D2
# add star to path
module add Alignment/STAR/2.7.0e
# B6 mm10
module add Alignment/STAR/2.7.0e; STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_B6mm10/ --genomeFastaFiles references/Mus_musculus.GRCm38.dna.toplevel.fa --outFileNamePrefix references/genome/star2.7.0e_B6mm10/B6mm10_ &
# D2
nohup STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_D2/ --genomeFastaFiles references/Mus_musculus_dba2j.DBA_2J_v1.dna.toplevel.fa --outFileNamePrefix references/genome/star2.7.0e_D2/D2_ > references/genome/star2.7.0e_D2/indexingD2.txt &
# B6 mm9
nohup STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_B6mm9/ --genomeFastaFiles references/Mus_musculus.NCBIM37.67.dna_rm.toplevel.fa --outFileNamePrefix references/genome/star2.7.0e_B6mm9/B6mm9_ > references/genome/star2.7.0e_B6mm9/indexingB6mm9.out 2> references/genome/star2.7.0e_B6mm9/indexingB6mm9.err &
# [waiting]

# Align with 0, 1, 2, 3, 4, or 5 mismatches

# path for output directory
##outdir=/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome
outdir=/home/ngobet/projects/BXD/data/TMP
### Problem with permission if output is in the NAS Franken (F)
# generate list of onliner commands
for ref in D2 B6_mm10 B6_mm9
do
	# all samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz; do
	# take only the parental and F1 samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*.fastq.gz; do
	# every sample except parental and F1
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/[^BD][^BD]*.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $fqfile
	##echo $sample
	# with 0 mismatch
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 99 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_0mismatch.out 2> $outdir/$ref/align\_$sample\_0mismatch.err " >> scripts/listcommands_alignRNAongenome_0mismatch.sh
	# with 1 mismatch maximum
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_1mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 98 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_1mismatch.out 2> $outdir/$ref/align\_$sample\_1mismatch.err " >> scripts/listcommands_alignRNAongenome_1mismatch.sh
	# with 2 mismatches maximum
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_2mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 97 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_2mismatch.out 2> $outdir/$ref/align\_$sample\_2mismatch.err " >> scripts/listcommands_alignRNAongenome_2mismatch.sh
	# with 3 mismatches maximum
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_3mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 96 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_3mismatch.out 2> $outdir/$ref/align\_$sample\_3mismatch.err " >> scripts/listcommands_alignRNAongenome_3mismatch.sh
	# with 4 mismatches maximum
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_4mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 95 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_4mismatch.out 2> $outdir/$ref/align\_$sample\_4mismatch.err " >> scripts/listcommands_alignRNAongenome_4mismatch.sh
	# with 5 mismatches maximum
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_5mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 94 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_5mismatch.out 2> $outdir/$ref/align\_$sample\_5mismatch.err " >> scripts/listcommands_alignRNAongenome_5mismatch.sh
	done
done

# Alignments options explanations:
# Decrease the score for insertions and deletions: --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 (so that alignment with insertions or deletions will have a AS at most of 99-10=89 so that I can filter it out)
# out filter insertions and deletions alignment --outFilterScoreMin 99 (score - (1*number of mismatch) 
# mismatch penalty is hardcoded at -1 (source: https://github.com/alexdobin/STAR/issues/140)

##nohup time bash listcommands.sh &
# launch the counting commands in parallel
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_0mismatch.sh 1> $outdir/0mismatch_parallel.out 2> $outdir/0mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_1mismatch.sh 1> $outdir/1mismatch_parallel.out 2> $outdir/1mismatch_parallel.err &

nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_2mismatch.sh 1> $outdir/2mismatch_parallel.out 2> $outdir/2mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_3mismatch.sh 1> $outdir/3mismatch_parallel.out 2> $outdir/3mismatch_parallel.err &

# testing how to combine nohup, parallel, and sleep commands
##sleep 0.002h; nohup parallel -j 1 <<< "echo 'Time is up'" 1> $outdir/test.out 2> $outdir/test.err &

nohup parallel -j 2 < scripts/listcommands_alignRNAongenome_4mismatch.sh 1> $outdir/4mismatch_parallel.out 2> $outdir/4mismatch_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_5mismatch.sh 1> $outdir/5mismatch_parallel.out 2> $outdir/5mismatch_parallel.err &

# copy files to F NAS
nohup mv /home/ngobet/projects/BXD/data/TMP/B6_mm9/*.bam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/B6_mm9 &
nohup mv /home/ngobet/projects/BXD/data/TMP/B6_mm10/* /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/B6_mm10 &
nohup mv /home/ngobet/projects/BXD/data/TMP/D2/* /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2 &

# index bam
/home/ngobet/software/samtools-1.9/bin/samtools index /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_0mismatchAligned.sortedByCoord.out.bam

# convert sam to bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_0mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_0mismatchAligned.sortedByCoord.out.bam

# test if correctly filtered (do we have some alignment with mismatche(s)

cut -f 10 data/transcriptome/2_mapping_STAR/OnGenome/D2/*1mismatch*sam | grep -c "N"
#50074
ngobet@franken-EESAA-server:/mnt/nas/BXD$ cat data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_1mismatchLog.final.out
#Number of input reads |       38032918
bc -l <<< "50074/38032918*100"
#.13165963232166409100
# Conclusion: there are reads with N gained but none of the uniquely mapped alignment has a real mismatch

/home/ngobet/software/samtools-1.9/bin/samtools view -o /home/ngobet/projects/BXD/data/TMP/D2/B61nsd_5mismatchAligned.sortedByCoord.out.sam /home/ngobet/projects/BXD/data/TMP/D2/B61nsd_5mismatchAligned.sortedByCoord.out.bam
head /home/ngobet/projects/BXD/data/TMP/D2/B61nsd_5mismatchAligned.sortedByCoord.out.sam
grep -v "nM:i:0" /home/ngobet/projects/BXD/data/TMP/D2/B61nsd_5mismatchAligned.sortedByCoord.out.sam
grep -vc "nM:i:1" /home/ngobet/projects/BXD/data/TMP/D2/B61nsd_5mismatchAligned.sortedByCoord.out.sam



# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on B6 mm10
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*0mismatchLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*0mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*1mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*2mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*3mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*4mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*5mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5

paste tmp_sample tmp_uniq* > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10_0mismatch.tsv
rm tmp_*

# retrieve statistics for alignment on D2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*0mismatchLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*0mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*1mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*2mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*3mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*4mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*5mismatchLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5

paste tmp_sample tmp_uniq* > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsD2_0mismatch.tsv
rm tmp_*
rm headers















# CHECK if really exactly matched

# convert sam to bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_0mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_0mismatchAligned.sortedByCoord.out.bam
# -h for include header

/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151518Aligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151518Aligned.sortedByCoord.out.bam
tail /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151518Aligned.sortedByCoord.out.sam

# testing command to have actually exact matches only!
nohup STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/05nsd\_0mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/05nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 > /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_0mismatch.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_0mismatch.err & 

# [copy]
##/home/ngobet/software/samtools-1.9/bin/samtools view -h -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/05nsd_0mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/05nsd_0mismatchAligned.sortedByCoord.out.bam
##tail /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/05nsd_0mismatchAligned.sortedByCoord.out.sam
##grep -v "101M" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/05nsd_0mismatchAligned.sortedByCoord.out.sam

/home/ngobet/software/samtools-1.9/bin/samtools view -h -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.bam
tail /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
cat /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatch*final*
#Total number of reads: 32904539
#Uniquely mapped reads number: 18323441
#Uniquely mapped reads % 55.69%
grep -c "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#18325171
grep -c "NH:i:1[^0]" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#18323441
grep "NH:i:1[^0]" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -c "101M"
#18179017
grep "NH:i:1[^0]" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep "101M" | grep -c "nM:i:0"
#18179017
grep "NH:i:1[^0]" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -c "nM:i:0"
#18323441
grep "NH:i:1[^0]" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep "nM:i:0" | grep -cv "101M"


grep "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -c "AS:i:99"
#18130772
grep "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -c "101M"
#18180747
grep "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -v "101M"
grep "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -c "nM:i:0"
#18325135
grep "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep "nM:i:0" | grep -c "HI:i:1"
#18323779
grep "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep "nM:i:0" | grep -v "HI:i:1"

##/home/ngobet/projects/BXD/data/TMP



# FIND THE CORRECT PARAMETERS

STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd\_0mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/DB1nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 > /home/ngobet/projects/BXD/data/TMP/D2/align\_DB1nsd\_0mismatch.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_DB1nsd\_0mismatch.err 

/home/ngobet/software/samtools-1.9/bin/samtools view -h -o /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.bam


grep -c "NH:i:1[^0]" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#18323441
grep "NH:i:1[^0]" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -v "101M"
grep "NH:i:1[^0]" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam | grep -c "101M"
#18179017
# check there are no multimappers
grep -v "NH:i:1[^0]" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#none
# check there are no insertion and deletion
grep -v "101M" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#there are some

grep -c "101M" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#18179017
grep -c "AS:i:99" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchAligned.sortedByCoord.out.sam
#18129082



# trying to change the penalty for insertions and deletions
outdir=/home/ngobet/projects/BXD/data/TMP
module add Alignment/STAR/2.7.0e
STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd\_0mismatchTEST --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/DB1nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > /home/ngobet/projects/BXD/data/TMP/D2/align\_DB1nsd\_0mismatchTEST.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_DB1nsd\_0mismatchTEST.err 

# added parameters
--scoreGap -10
--scoreDelOpen -10
--scoreInsOpen -10

#convert to sam
/home/ngobet/software/samtools-1.9/bin/samtools view -h -o /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchTESTAligned.sortedByCoord.out.sam /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchTESTAligned.sortedByCoord.out.bam

grep -v "101M" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchTESTAligned.sortedByCoord.out.sam
#there are some

# testing with only taking reads with alignment score 99
--outFilterScoreMin 99
STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd\_0mismatchTEST --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/DB1nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 99 > /home/ngobet/projects/BXD/data/TMP/D2/align\_DB1nsd\_0mismatchTEST.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_DB1nsd\_0mismatchTEST.err 

#convert to sam
/home/ngobet/software/samtools-1.9/bin/samtools view -h -o /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchTESTAligned.sortedByCoord.out.sam /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchTESTAligned.sortedByCoord.out.bam

grep -v "101M" /home/ngobet/projects/BXD/data/TMP/D2/DB1nsd_0mismatchTESTAligned.sortedByCoord.out.sam
#No alignment has insertions and deletions





--alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMatchNminOverLread 1.0
## source: https://groups.google.com/forum/#!searchin/rna-star/insertions%7Csort:date/rna-star/R-CBNFeLlgU/lgeOBIkfCAAJ
--alignIntronMax 1
## source: https://groups.google.com/forum/#!topic/rna-star/a5cZf7CateQ




# Filtering output parameters
# 19 June 2019

--outFilterMultimapNmax 1
# default: 10 int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as ”mapped to too many loci” in the Log.final.out .

--outFilterMismatchNmax 0
# default: 10 int: alignment will be output only if it has no more mismatches than this value.

--outFilterMatchNminOverLread 1.0
# default: 0.66 real: same as outFilterMatchNmin, but normalized to the read length (sum of mates’ lengths for paired-end reads).	--outFilterMatchNmin default: 0 int: alignment will be output only if the number of matched bases is higher than or equal to this value.

--alignIntronMax 1
# default: 0 maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins

--alignEndsType EndToEnd
# default: Local string: type of read ends alignment
	# Local standard local alignment with soft-clipping allowed
	# EndToEnd force end-to-end read alignment, do not soft-clip
	# Extend5pOfRead1 fully extend only the 5p of the read1, all other ends: local alignment
	# Extend5pOfReads12 fully extend only the 5p of the both read1 and read2, all other ends:local alignment



# not important
--outFilterMismatchNoverLmax 0
--outFilterMismatchNoverReadLmax 1.0
--outFilterIntronMotifs
outFilterIntronStrands


# testing outFilterMatchNmin parameters
# 1.0
module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/05nsd\_3mismatchTEST1 --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/05nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 96 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_3mismatchTEST1.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_3mismatchTEST1.err 
# 0.97
module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/05nsd\_3mismatchTEST097 --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/05nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 0.97 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 96 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_3mismatchTEST097.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_3mismatchTEST097.err 
# 0.96
module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_D2 --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/D2/05nsd\_3mismatchTEST096 --runThreadN 8 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/05nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 0.96 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 96 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_3mismatchTEST096.out 2> /home/ngobet/projects/BXD/data/TMP/D2/align\_05nsd\_3mismatchTEST096.err 


# move the files
mv /home/ngobet/projects/BXD/data/TMP/D2/*TEST* /mnt/nas/BXD/

# index bam
/home/ngobet/software/samtools-1.9/bin/samtools index /mnt/nas/BXD/*TEST1*.bam
/home/ngobet/software/samtools-1.9/bin/samtools index /mnt/nas/BXD/*TEST097*.bam
/home/ngobet/software/samtools-1.9/bin/samtools index /mnt/nas/BXD/*TEST096*.bam

# convert sam to bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/05nsd_3mismatchTEST096Aligned.sortedByCoord.out.sam /mnt/nas/BXD/05nsd_3mismatchTEST096Aligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/05nsd_3mismatchTEST1Aligned.sortedByCoord.out.sam /mnt/nas/BXD/05nsd_3mismatchTEST1Aligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/05nsd_3mismatchTEST097Aligned.sortedByCoord.out.sam /mnt/nas/BXD/05nsd_3mismatchTEST097Aligned.sortedByCoord.out.bam

grep "nM:i:3" /mnt/nas/BXD/*.sam

/home/ngobet/software/bedtools2/bin/intersectBed -v -bed -a 05nsd_3mismatchTEST096Aligned.sortedByCoord.out.bam -b 05nsd_3mismatchTEST097Aligned.sortedByCoord.out.bam > TMPdiffAlignment096_097



# 27 June 2019

# testing parameters


# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP/TEST
# generate list of onliner commands
for ref in D2 B6_mm10
do
	# on B6 and D2 samples, 1st replicate, in liver (smaller bam)
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L[BD][B6]1nsd.fastq.gz; do
	# on B6, 1st replicate, in liver (smaller bam)
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/LB61nsd.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $fqfile
	##echo $sample
	# default (remove --outFilterMismatchNmax 0 --alignIntronMax 1 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 --alignEndsType EndToEnd --outFilterMultimapNmax 1)
	##spe="free"
	##echo $spe
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# as Ioannis wants
	spe="exactunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	spe="1mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	spe="2mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	spe="3mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	spe="4mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	spe="5mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# control mismatch
	##spe="0mismatch"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMismatchNmax 0 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##spe="1mismatch"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMismatchNmax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##spe="2mismatch"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMismatchNmax 2 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##spe="3mismatch"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMismatchNmax 3 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# Indels control
	# no insertions and deletions
	##spe="noindels"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# no introns
	##spe="nointrons1"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --alignIntronMax 1 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##spe="nointrons2"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --alignIntronMax 1 --scoreGap -40 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##spe="nointrons3"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreGap -40 > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# no softclip
	##spe="nosoftclip"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_TEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --alignEndsType EndToEnd > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# no all
	spe="onlyM"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0mismatchTEST$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd > $outdir/$ref/align\_$sample\_TEST$spe.out 2> $outdir/$ref/align\_$sample\_TEST$spe.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# with 0 mismatch
	#with or without score restriction (--outFilterScoreMin)
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0mismatchTESTnoAS --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_0mismatchTESTnoAS.out 2> $outdir/$ref/align\_$sample\_0mismatchTESTnoAS.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0mismatchTESTAS99 --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 99 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_0mismatchTESTAS99.out 2> $outdir/$ref/align\_$sample\_0mismatchTESTAS99.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0mismatchTESTAS98 --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 98 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_0mismatchTESTAS98.out 2> $outdir/$ref/align\_$sample\_0mismatchTESTAS98.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# remove --outFilterMatchNminOverLread 1.0
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_0mismatchTESTnoASnooutMatchOver --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --alignIntronMax 1 --outFilterMultimapNmax 1 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_0mismatchTESTnoASnooutMatchOver.out 2> $outdir/$ref/align\_$sample\_0mismatchTESTnoASnooutMatchOver.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# with 1 mismatch maximum
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_1mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 98 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_1mismatch.out 2> $outdir/$ref/align\_$sample\_1mismatch.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# with 2 mismatches maximum
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_2mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 97 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_2mismatch.out 2> $outdir/$ref/align\_$sample\_2mismatch.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# with 3 mismatches maximum
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_3mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 96 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_3mismatch.out 2> $outdir/$ref/align\_$sample\_3mismatch.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# with 4 mismatches maximum
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_4mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 95 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_4mismatch.out 2> $outdir/$ref/align\_$sample\_4mismatch.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	# with 5 mismatches maximum
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_5mismatch --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMatchNminOverLread 1.0 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterScoreMin 94 --scoreGap -10 --scoreDelOpen -10 --scoreInsOpen -10 > $outdir/$ref/align\_$sample\_5mismatch.out 2> $outdir/$ref/align\_$sample\_5mismatch.err " >> scripts/listcommands_alignRNAongenome_TEST.sh
	done
done

# running test alignments
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_TEST.sh 1> $outdir/TEST_parallel.out 2> $outdir/TEST_parallel.err &

# move files
###remove the incorrect 0mismatch before to move!
mv /home/ngobet/projects/BXD/data/TMP/TEST/D2/*TEST* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2 ####
mv /home/ngobet/projects/BXD/data/TMP/TEST/B6_mm10/*TEST* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6_mm10

# index bam
ls /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/*.bam | xargs -n1 -P1 /home/ngobet/software/samtools-1.9/bin/samtools index

/home/ngobet/software/samtools-1.9/bin/samtools index /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.bam
# help from: https://www.biostars.org/p/170522/

# convert sam to bam (/home/ngobet/software/samtools-1.9/bin/samtools view -o)
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTAS98Aligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTAS98Aligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTnoASAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTnoASAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTAS99Aligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTAS99Aligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LDB1nsd_0mismatchTESTfreeAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LDB1nsd_0mismatchTESTfreeAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTfreeAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTfreeAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTnoindelsAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTnoindelsAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST0mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST0mismatchAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST1mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST1mismatchAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST2mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST2mismatchAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST3mismatchAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST3mismatchAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTnosoftclipAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTnosoftclipAligned.sortedByCoord.out.bam

/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.bam


# checking filtering

# Do they have no insertion and deletion?
grep -cv "101M" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/*.sam
#0 for AS99 and AS98, 163171 for noAS


# Do they have the required alignment score quality?
grep -cv "AS:i:99" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/*.sam
#AS98:13032, AS99:0, noAS:176203 => works

# 
grep "Average" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/*final*


# Drop in quality due to N in sequences?
grep -v "AS:i:99" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTAS98Aligned.sortedByCoord.out.sam | cut -f 10 | grep -c N
#7144/13032
grep -v "AS:i:99" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTnoASAligned.sortedByCoord.out.sam | cut -f 10 | grep -c N
#7197/176203

awk ' $10 { print length }' /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_0mismatchTESTnoASAligned.sortedByCoord.out.sam 



# Check rigorously filtering

# No mismatch? grep -cv "nM:i:0" should be 0
grep -cv "nM:i:0" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST0mismatchAligned.sortedByCoord.out.sam
#0. This works.
# No more than 1 mismatch? grep -cv "nM:i:[01]" should be 0
grep -cv "nM:i:[01]" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST1mismatchAligned.sortedByCoord.out.sam
#0
# No more than 2 mismatches? grep -c "nM:i:[012]" should be 0
grep -cv "nM:i:[012]" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST2mismatchAligned.sortedByCoord.out.sam
#0
# No more than 3 mismatches? grep -c "nM:i:[0123]" should be 0
grep "nM:i:3" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST3mismatchAligned.sortedByCoord.out.sam | head
# yes there are some alignments with 3 matches
grep -cv "nM:i:[0123]" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TEST3mismatchAligned.sortedByCoord.out.sam
#0 it works
# Do they have no insertion and deletion? grep -c "[DI]" on the cigar string should be 0 and insertions and deletions rates in *final* are 0.
cut -f 6 /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTnoindelsAligned.sortedByCoord.out.sam | grep -c "[DI]"
#0. It worked.
# Do they have no introns? *SJ*.tab should be 0 (empty) and grep "Number of splices: Total" in *final* should be 0
#This is ok for the 3 commands to remove introns. => keeping --alignIntronMax 1 (Proposed by the programmer)
# Do they have no softclipping? ### would it appear in the CIGAR?
cut -f 6 /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTnosoftclipAligned.sortedByCoord.out.sam | grep -c "S"
#0 It works.


# ultime test
# exact?
# no mismatch
grep -cv "nM:i:0" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.sam
#0
# no indels or softclip
grep -cv "101M" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.sam
#0
# no introns
grep "Number of splices: Total" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactunique*final*
#0
# uniq?
grep -cv "NH:i:1" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.sam
#0
# They are exact and uniq. :)


# Out of curiosity, how many have alignment score (AS) less than 99?
grep -cv "AS:i:99" /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2/LB61nsd_TESTexactuniqueAligned.sortedByCoord.out.sam
#72070



# retrieve a few values from the organized test
# get the name of test
grep "Uniquely mapped reads %" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/TEST/D2/old/LB61nsd*final* | cut -d _ -f 4 | cut -d L -f 1 > tmp_name
# get the uniquely mapped percentage
grep "Uniquely mapped reads %" /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/TEST/D2/old/LB61nsd*final* | cut -f 2 | cut -d % -f 1 > tmp_uniqpercentage

paste tmp_name tmp_uniqpercentage > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsTEST_LB61nsd_onD2.tsv
rm tmp_*


####################################################################################################################


# alignments

# set working directory
cd /mnt/nas/BXD

# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP
# generate list of onliner commands
##for ref in D2 B6mm10
for ref in B6mm10_primaryassembly
do
	# for all samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz; do
	# take only the parental and F1 samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*.fastq.gz; do
	# for BXD samples (every sample except parental and F1)
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/[^BD][^BD]*.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $fqfile
	##echo $sample
	spe="exactunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="1mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="2mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="3mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="4mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="5mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##spe="6mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 6 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="10mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 10 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	done
done

# running parental and F1 alignments
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_exactunique.sh 1> $outdir/exactunique_parallel.out 2> $outdir/exactunique_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_1mismatchunique.sh 1> $outdir/1mismatchunique_parallel.out 2> $outdir/1mismatchunique_parallel.err &

nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_2mismatchunique.sh 1> $outdir/2mismatchunique_parallel.out 2> $outdir/2mismatchunique_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_3mismatchunique.sh 1> $outdir/3mismatchunique_parallel.out 2> $outdir/3mismatchunique_parallel.err &

nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_4mismatchunique.sh 1> $outdir/4mismatchunique_parallel.out 2> $outdir/4mismatchunique_parallel.err &
nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_5mismatchunique.sh 1> $outdir/5mismatchunique_parallel.out 2> $outdir/5mismatchunique_parallel.err &

nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_10mismatchunique.sh 1> $outdir/10mismatchunique_parallel.out 2> $outdir/10mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_6mismatchunique.sh 1> $outdir/6mismatchunique_parallel.out 2> $outdir/6mismatchunique_parallel.err &


# separating NSD and SD
##split -l 2 -d -a 3 scripts/listcommandsBXD_alignRNAongenome_exactunique.sh
##cat x*[02468] > scripts/listcommandsBXD_alignRNAongenome_exactuniqueNSD.sh
##cat x*[13579] > scripts/listcommandsBXD_alignRNAongenome_exactuniqueSD.sh
##rm x*
##split -l 2 -d -a 3 scripts/listcommandsBXD_alignRNAongenome_1mismatchunique.sh
##cat x*[02468] > scripts/listcommandsBXD_alignRNAongenome_1mismatchuniqueNSD.sh
##cat x*[13579] > scripts/listcommandsBXD_alignRNAongenome_1mismatchuniqueSD.sh
##rm x*
##split -l 2 -d -a 3 scripts/listcommandsBXD_alignRNAongenome_2mismatchunique.sh
##cat x*[02468] > scripts/listcommandsBXD_alignRNAongenome_2mismatchuniqueNSD.sh
##cat x*[13579] > scripts/listcommandsBXD_alignRNAongenome_2mismatchuniqueSD.sh
##rm x*
##split -l 2 -d -a 3 scripts/listcommandsBXD_alignRNAongenome_3mismatchunique.sh
##cat x*[02468] > scripts/listcommandsBXD_alignRNAongenome_3mismatchuniqueNSD.sh
##cat x*[13579] > scripts/listcommandsBXD_alignRNAongenome_3mismatchuniqueSD.sh
##rm x*
##split -l 2 -d -a 3 scripts/listcommandsBXD_alignRNAongenome_4mismatchunique.sh
##cat x*[02468] > scripts/listcommandsBXD_alignRNAongenome_4mismatchuniqueNSD.sh
##cat x*[13579] > scripts/listcommandsBXD_alignRNAongenome_4mismatchuniqueSD.sh
##rm x*
##split -l 2 -d -a 3 scripts/listcommandsBXD_alignRNAongenome_5mismatchunique.sh
##cat x*[02468] > scripts/listcommandsBXD_alignRNAongenome_5mismatchuniqueNSD.sh
##cat x*[13579] > scripts/listcommandsBXD_alignRNAongenome_5mismatchuniqueSD.sh
##rm x*

# running BXD alignments
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_exactunique.sh 1> $outdir/exactunique_parallel.out 2> $outdir/exactunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_1mismatchunique.sh 1> $outdir/1mismatchunique_parallel.out 2> $outdir/1mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_2mismatchunique.sh 1> $outdir/2mismatchunique_parallel.out 2> $outdir/2mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_3mismatchunique.sh 1> $outdir/3mismatchunique_parallel.out 2> $outdir/3mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_4mismatchunique.sh 1> $outdir/4mismatchunique_parallel.out 2> $outdir/4mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_5mismatchunique.sh 1> $outdir/5mismatchunique_parallel.out 2> $outdir/5mismatchunique_parallel.err &
# NSD BXD samples
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_exactuniqueNSD.sh 1> /home/ngobet/projects/BXD/data/TMP/exactuniqueNSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/exactuniqueNSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_1mismatchuniqueNSD.sh 1> /home/ngobet/projects/BXD/data/TMP/1mismatchuniqueNSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/1mismatchuniqueNSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_2mismatchuniqueNSD.sh 1> /home/ngobet/projects/BXD/data/TMP/2mismatchuniqueNSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/2mismatchuniqueNSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_3mismatchuniqueNSD.sh 1> /home/ngobet/projects/BXD/data/TMP/3mismatchuniqueNSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/3mismatchuniqueNSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_4mismatchuniqueNSD.sh 1> /home/ngobet/projects/BXD/data/TMP/4mismatchuniqueNSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/4mismatchuniqueNSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_5mismatchuniqueNSD.sh 1> /home/ngobet/projects/BXD/data/TMP/5mismatchuniqueNSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/5mismatchuniqueNSD_parallel.err
# SD BXD samples
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_exactuniqueSD.sh 1> /home/ngobet/projects/BXD/data/TMP/exactuniqueSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/exactuniqueSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_1mismatchuniqueSD.sh 1> /home/ngobet/projects/BXD/data/TMP/1mismatchuniqueSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/1mismatchuniqueSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_2mismatchuniqueSD.sh 1> /home/ngobet/projects/BXD/data/TMP/2mismatchuniqueSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/2mismatchuniqueSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_3mismatchuniqueSD.sh 1> /home/ngobet/projects/BXD/data/TMP/3mismatchuniqueSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/3mismatchuniqueSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_4mismatchuniqueSD.sh 1> /home/ngobet/projects/BXD/data/TMP/4mismatchuniqueSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/4mismatchuniqueSD_parallel.err
##nohup parallel -j 1 < scripts/listcommandsBXD_alignRNAongenome_5mismatchuniqueSD.sh 1> /home/ngobet/projects/BXD/data/TMP/5mismatchuniqueSD_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/5mismatchuniqueSD_parallel.err

# Retrieve info

# parental and F1
# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
echo "Uniq6" >> tmp_headers
echo "Uniq10" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on D2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*6mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq6
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*[BD]*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq6 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsD2.tsv
rm tmp_*

# retrieve statistics for alignment on B6 mm10
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*6mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq6
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6_mm10/*[BD]*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq6 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10_B6D2F1.tsv
rm tmp_*

rm headers

# all samples
# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
echo "Uniq10" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on D2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsD2.tsv
rm tmp_*

# retrieve statistics for alignment on B6 mm10
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10/*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10.tsv
rm tmp_*

rm headers

################################


# check percentage of reads with many mismatches in alignment to transcriptome reference

# index bam
##/home/ngobet/software/samtools-1.9/bin/samtools index /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_0mismatchAligned.sortedByCoord.out.bam

# convert sam to bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151541Aligned.sortedByCoord.out.sam data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151541Aligned.sortedByCoord.out.bam

grep "nM:i:10" data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151541Aligned.sortedByCoord.out.sam
#171064
grep -vc "nM:i:[01]" data/transcriptome/2_mapping_STAR/2-pass/D2/GSM3151541Aligned.sortedByCoord.out.sam
#5124789

#(5124789+171064)/40819200*100 = 12.97393
# 13% of reads contain more than 1 mismatch.

################################

# test influence of leaving all patches from the genome or only main chromosomes

# find files without all the patches
ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gzip -d *.gz
# make a file without patches
sed '/JH584299.1 dna:scaffold scaffold:GRCm38:JH584299.1:1:953012:1 REF/q' references/Mus_musculus.GRCm38.dna.primary_assembly.fa > references/Mus_musculus.GRCm38.dna.primary_assembly_noscaffolds.faTMP
head -n -1 references/Mus_musculus.GRCm38.dna.primary_assembly_noscaffolds.faTMP > references/Mus_musculus.GRCm38.dna.primary_assembly_noscaffolds.fa
rm references/Mus_musculus.GRCm38.dna.primary_assembly_noscaffolds.faTMP
# help from: https://unix.stackexchange.com/questions/170661/awk-print-lines-from-the-first-line-until-match-word
# help from: https://stackoverflow.com/questions/4881930/remove-the-last-line-from-a-file-in-bash
# make a file without patches and Y and MT chromosomes
grep ">" -n references/Mus_musculus.GRCm38.dna.primary_assembly.fa
#1:>1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
#38969201:>9 dna:chromosome chromosome:GRCm38:9:1:124595110:1 REF
#41045788:>MT dna:chromosome chromosome:GRCm38:MT:1:16299:1 REF
#41046061:>X dna:chromosome chromosome:GRCm38:X:1:171031299:1 REF
#43896584:>Y dna:chromosome chromosome:GRCm38:Y:1:91744698:1 REF
#45425664:>JH584299.1 dna:scaffold scaffold:GRCm38:JH584299.1:1:953012:1 REF
#45514596:>JH584295.1 dna:scaffold scaffold:GRCm38:JH584295.1:1:1976:1 REF
sed -n '1,41045787p;41046061,43896583p;43896584q' references/Mus_musculus.GRCm38.dna.primary_assembly.fa > references/Mus_musculus.GRCm38.dna.major_chromosomes.fa
# help from: https://stackoverflow.com/questions/83329/how-can-i-extract-a-predetermined-range-of-lines-from-a-text-file-on-unix

# generate genome + transcriptome
# Generate the reference genomes (without annotation) for star

# create directories
mkdir -p references/genome/star2.7.0e_B6mm10_primaryassembly
mkdir -p references/genome/star2.7.0e_B6mm10_primaryassembly_noscaffolds
mkdir -p references/genome/star2.7.0e_B6mm10_majorchromosomes

# add star to path
module add Alignment/STAR/2.7.0e

# B6 mm10 primary assembly (without patches but with scaffolds, recommended in STAR manual)
nohup STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_B6mm10_primaryassembly/ --genomeFastaFiles references/Mus_musculus.GRCm38.dna.primary_assembly.fa --outFileNamePrefix references/genome/star2.7.0e_B6mm10_primaryassembly/B6mm10_primaryassembly_ > references/genome/star2.7.0e_B6mm10_primaryassembly/indexingB6mm10_primaryassembly.out 2> references/genome/star2.7.0e_B6mm10_primaryassembly/indexingB6mm10_primaryassembly.err &
# B6 mm10 primary assembly without patches
nohup STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_B6mm10_primaryassembly_noscaffolds/ --genomeFastaFiles references/Mus_musculus.GRCm38.dna.primary_assembly_noscaffolds.fa --outFileNamePrefix references/genome/star2.7.0e_B6mm10_primaryassembly_noscaffolds/B6mm10_primaryassembly_noscaffolds_ > references/genome/star2.7.0e_B6mm10_primaryassembly_noscaffolds/indexingB6mm10_primaryassembly_noscaffolds.out 2> references/genome/star2.7.0e_B6mm10_primaryassembly_noscaffolds/indexingB6mm10_primaryassembly_noscaffolds.err &
# B6 mm10 major chromosomes
nohup STAR --runThreadN 4 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_B6mm10_majorchromosomes/ --genomeFastaFiles references/Mus_musculus.GRCm38.dna.major_chromosomes.fa --outFileNamePrefix references/genome/star2.7.0e_B6mm10_majorchromosomes/B6mm10_majorchromosomes_ > references/genome/star2.7.0e_B6mm10_majorchromosomes/indexingB6mm10_majorchromosomes.out 2> references/genome/star2.7.0e_B6mm10_majorchromosomes/indexingB6mm10_majorchromosomes.err &


# test on 1 or few samples
 
# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP

# generate list of oneliner commands
# take only the parental nsd samples
for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz; do
sample=$(basename -s .fastq.gz $fqfile)
##echo $fqfile
##echo $sample
	for ref in B6mm10_majorchromosomes B6mm10_primaryassembly B6mm10_primaryassembly_noscaffolds
	do
	spe="exactunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="1mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="2mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="3mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="4mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="5mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="10mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 10 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	done
done

# running BXD alignments
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_exactunique.sh 1> $outdir/exactunique_parallel.out 2> $outdir/exactunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_1mismatchunique.sh 1> $outdir/1mismatchunique_parallel.out 2> $outdir/1mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_2mismatchunique.sh 1> $outdir/2mismatchunique_parallel.out 2> $outdir/2mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_3mismatchunique.sh 1> $outdir/3mismatchunique_parallel.out 2> $outdir/3mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_4mismatchunique.sh 1> $outdir/4mismatchunique_parallel.out 2> $outdir/4mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_5mismatchunique.sh 1> $outdir/5mismatchunique_parallel.out 2> $outdir/5mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenomeTEST_10mismatchunique.sh 1> $outdir/5mismatchunique_parallel.out 2> $outdir/10mismatchunique_parallel.err & 

# length D2 genome with scaffolds
awk '{sum+=$2} END {print sum}' references/genome/star2.7.0e_D2/chrNameLength.txt
#2.60615e+09

# length D2 genome without scaffolds (chromosomes 1 to 19 and X)
head -n 20 references/genome/star2.7.0e_D2/chrNameLength.txt | awk '{sum+=$2} END {print sum}'
#2.57879e+09


# length D2 scaffolds
bc <<< "2606145858-2578793109"
#27352749

# get length of B6 mm10 genome (with/without patches/scaffolds/chromosomes Y MT)
awk '{sum+=$2} END {print sum}' references/genome/star2.7.0e_B6_mm10/chrNameLength.txt
#1.25513e+10
awk '{sum+=$2} END {print sum}' references/genome/star2.7.0e_B6mm10_primaryassembly/chrNameLength.txt
#2.73087e+09
awk '{sum+=$2} END {print sum}' references/genome/star2.7.0e_B6mm10_primaryassembly_noscaffolds/chrNameLength.txt
#2.72554e+09
awk '{sum+=$2} END {print sum}' references/genome/star2.7.0e_B6mm10_majorchromosomes/chrNameLength.txt
#2.63378e+09 (chromosomes 1 to 19 and X)


# Retrieve info
# parental and F1
# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
echo "Uniq10" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on B6 mm10 primary assembly
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 7 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*1mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*2mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*3mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*4mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*5mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly/*[BD]*[12]nsd*10mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10primaryassembly_B6D2.tsv
rm tmp_*

# retrieve statistics for alignment on B6 mm10 primary assembly without scaffolds
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 7 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*1mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*2mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*3mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*4mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*5mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_primaryassembly_noscaffolds/*[BD]*[12]nsd*10mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10primaryassemblynoscaffolds_B6D2.tsv
rm tmp_*

# retrieve statistics for alignment on B6 mm10 majorchromosomes
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 7 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*1mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*2mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*3mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*4mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*5mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes/*[BD]*[12]nsd*10mismatchuniqueLog.final.out | cut -d / -f 7 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10majorchromosomes_B6D2.tsv
rm tmp_*


rm headers


################################

# test influence of leaving scaffolds or only main chromosomes for D2
grep ">" -n references/Mus_musculus_dba2j.DBA_2J_v1.dna.toplevel.fa | head -n 30
#1:>1 dna:chromosome chromosome:DBA_2J_v1:1:1:192928109:1 REF
#39233561:>19 dna:chromosome chromosome:DBA_2J_v1:19:1:58689003:1 REF
#40211713:>X dna:chromosome chromosome:DBA_2J_v1:X:1:166092115:1 REF
#42979916:>KV417258.1 dna:scaffold scaffold:DBA_2J_v1:KV417258.1:1:359851:1 REF
sed -n '1,42979915p;42979916q' references/Mus_musculus_dba2j.DBA_2J_v1.dna.toplevel.fa > references/Mus_musculus_dba2j.DBA_2J_v1.dna.major_chromosomes.fa
# test influence of leaving scaffolds but not Y and MT chromosomes for B6 mm10
sed -n '1,41045787p;41046061,43896583p;45425664,45514629p' references/Mus_musculus.GRCm38.dna.primary_assembly.fa > references/Mus_musculus.GRCm38.dna.major_chromosomes_scaffolds.fa

# create directories
mkdir -p references/genome/star2.7.0e_D2_majorchromosomes
mkdir -p data/transcriptome/2_Mapping_STAR/OnGenome/TEST/D2_majorchromosomes
mkdir -p /home/ngobet/projects/BXD/data/TMP/D2_majorchromosomes

mkdir -p references/genome/star2.7.0e_B6mm10_majorchromosomes_scaffolds
mkdir -p data/transcriptome/2_Mapping_STAR/OnGenome/TEST/B6mm10_majorchromosomes_scaffolds
mkdir -p /home/ngobet/projects/BXD/data/TMP/B6mm10_majorchromosomes_scaffolds



# add star to path
module add Alignment/STAR/2.7.0e

# generate index genomes
nohup STAR --runThreadN 8 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_D2_majorchromosomes/ --genomeFastaFiles references/Mus_musculus_dba2j.DBA_2J_v1.dna.major_chromosomes.fa --outFileNamePrefix references/genome/star2.7.0e_D2_majorchromosomes/D2_majorchromosomes_ > references/genome/star2.7.0e_D2_majorchromosomes/indexingD2_majorchromosome.out 2> references/genome/star2.7.0e_D2_majorchromosomes/indexingD2_majorchromosome.err &

 
# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP

# generate list of oneliner commands
# take only the parental nsd samples
for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz; do
sample=$(basename -s .fastq.gz $fqfile)
##echo $fqfile
##echo $sample
	##for ref in D2_majorchromosomes B6mm10_majorchromosomes_scaffolds
	for ref in B6mm10_majorchromosomes_scaffolds
	do
	spe="exactunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="1mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="2mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="3mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="4mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="5mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	spe="10mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 10 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/TEST/$ref" >> scripts/listcommands_alignRNAongenomeTEST\_$spe.sh
	done
done

####################################################################################################################

# 19 July 2019

# alignments

# set working directory
cd /mnt/nas/BXD

# path for output directory
outdir=/home/ngobet/projects/BXD/data/TMP

# generate list of onliner commands
for ref in D2_majorchromosomes B6mm10_majorchromosomes
##for ref in B6mm10_majorchromosomes
do
	# for all samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz; do
	# take only the parental and F1 samples
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*.fastq.gz; do
	# take only the parental nsd samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz; do
	# for BXD samples (every sample except parental and F1)
	##for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/[^BD][^BD]*.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $fqfile
	##echo $sample
	##spe="exactunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##spe="1mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="2mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="3mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="4mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 4 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="5mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##spe="6mismatchunique"
	##echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 6 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	##echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	spe="10mismatchunique"
	echo "module add Alignment/STAR/2.7.0e; STAR --genomeDir references/genome/star2.7.0e_$ref --outFileNamePrefix $outdir/$ref/$sample\_$spe --runThreadN 8 --readFilesCommand zcat --readFilesIn $fqfile --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 10 --outFilterMultimapNmax 1 > $outdir/$ref/align\_$sample\_$spe.out 2> $outdir/$ref/align\_$sample\_$spe.err " >> scripts/listcommands_alignRNAongenome\_$spe.sh
	echo "mv /home/ngobet/projects/BXD/data/TMP/$ref/*$spe* /mnt/nas/BXD/data/transcriptome/2_Mapping_STAR/OnGenome/$ref" >> scripts/listcommands_alignRNAongenome\_$spe.sh
	done
done

# running parental and F1 alignments
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_exactunique.sh 1> /home/ngobet/projects/BXD/data/TMP/exactunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/exactunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_1mismatchunique.sh 1> /home/ngobet/projects/BXD/data/TMP/1mismatchunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/1mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_2mismatchunique.sh 1> /home/ngobet/projects/BXD/data/TMP/2mismatchunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/2mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_3mismatchunique.sh 1> /home/ngobet/projects/BXD/data/TMP/3mismatchunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/3mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_4mismatchunique.sh 1> /home/ngobet/projects/BXD/data/TMP/4mismatchunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/4mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_5mismatchunique.sh 1> /home/ngobet/projects/BXD/data/TMP/5mismatchunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/5mismatchunique_parallel.err &
##nohup parallel -j 1 < scripts/listcommands_alignRNAongenome_10mismatchunique.sh 1> /home/ngobet/projects/BXD/data/TMP/10mismatchunique_parallel.out 2> /home/ngobet/projects/BXD/data/TMP/10mismatchunique_parallel.err &

# Retrieve info
# parental and F1
# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
echo "Uniq10" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on B6 mm10 majorchromosomes
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_majorchromosomes/*[BD]*[12]nsd*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10majorchromosomes_B6D2.tsv
rm tmp_*

# retrieve statistics for alignment on D2 majorchromosomes
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2_majorchromosomes/*[BD]*[12]nsd*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsD2majorchromosomes_B6D2.tsv
rm tmp_*

rm headers


################################
# bed coverage 

# convert bam to bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_exactuniqueAligned.sortedByCoord.out.bam -cigar > B61nsd_onD2_0.bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnGenome/B6_mm10/B61nsd_exactuniqueAligned.sortedByCoord.out.bam -cigar > B61nsd_onB6mm10_0.bed

/home/ngobet/software/bedtools2/bin/coverageBed -a B61nsd_onD2_0.bed -b B61nsd_onB6mm10_0.bed > tmp.cov
/home/ngobet/software/bedtools2/bin/coverageBed -b B61nsd_onD2_0.bed -a B61nsd_onB6mm10_0.bed > tmpBA.cov

##
#20506164 B61nsd_onB6mm10_0.bed
# 17555054 B61nsd_onD2_0.bed
# 20506164 tmpBA.cov
# 17555054 tmp.cov
# 76122436 total
#Default Output:
#        After each entry in A, reports:
#          1) The number of features in B that overlapped the A interval.
#          2) The number of bases in A that had non-zero coverage.
#          3) The length of the entry in A.
#          4) The fraction of bases in A that had non-zero coverage.

# find region of interest
grep "0\.[^0]" tmp.cov

##############################

# 29 August 2019

# difference dna dna_sm and dna_rm

# hard vs soft-masked in STAR
https://www.biostars.org/p/290455/
https://groups.google.com/forum/#!topic/rna-star/2wdHXaPv_vU

# Which repeats are soft-masked?
https://www.biostars.org/p/152726/


ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus_dba2j/dna/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa.gz

# 20 August 2019

# add star to path
module add Alignment/STAR/2.7.0e

# generate index genomes
nohup STAR --runThreadN 8 --limitGenomeGenerateRAM 33524399488 --runMode genomeGenerate --genomeDir references/genome/star2.7.0e_B6mm10_primaryassembly/ --genomeFastaFiles references/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa --outFileNamePrefix references/genome/star2.7.0e_B6mm10_primaryassembly/B6mm10_primaryassembly_ > references/genome/star2.7.0e_B6mm10_primaryassembly/indexingB6mm10_primaryassembly.out 2> references/genome/star2.7.0e_B6mm10_primaryassembly/indexingB6mm10_primaryassembly.err &

############################

#PART 3: Counting reads per gene

#  2 September 2019
# required: htseq and samtools
# documentation htseq: https://htseq.readthedocs.io/en/release_0.10.0/counting.html

# set working directory:
cd /mnt/nas/BXD/

#input directory
indir=data/transcriptome/2_mapping_STAR/OnGenome

#output directory
outdir=data/transcriptome/3_counting_htseq/OnGenome
mkdir -p $outdir/B6mm10_primaryassembly/
mkdir -p $outdir/D2/


# generate list of commands to count reads per gene for 3 different references
for ref in B6mm10_primaryassembly D2
do
##echo $ref
	for bamfile in $indir/$ref/*_exactuniqueAligned.sortedByCoord.out.bam
	do
	sample=$(basename -s _exactuniqueAligned.sortedByCoord.out.bam $bamfile)
	##echo $sample
	if [ "$ref" == "D2" ]; then
	reffile=references/Mus_musculus_dba2j.DBA_2J_v1.94.gtf
	elif [ "$ref" == "B6mm10_primaryassembly" ]; then
	reffile=references/Mus_musculus.GRCm38.94.gtf
	fi
	echo "/home/ngobet/software/samtools-1.9/bin/samtools view -h $bamfile | htseq-count -q -s reverse -t exon -m union - $reffile 2>$outdir/$ref/$sample\_htseqcount.err 1>$outdir/$ref/$sample.count" >> /mnt/nas/BXD/$outdir/Counting_listcommands.sh
	done
done

# launch the counting commands in parallel
nohup parallel -j 4 < $outdir/Counting_listcommands.sh > $outdir/parallel.out 2> $outdir/parallel.err &

# group count file to have one file with gene names as rows and BXD line as columns.
###bsub -o results/transcriptome/3_counting_htseq/bsub_summary_counts.out -e results/transcriptome/3_counting_htseq/bsub_summary_counts.err -J bsub_summary_counts scripts/summary_counts.sh

# 3 September 2019
# merging individual count file per sample into one file per reference
bash scripts/MergeCount.sh

# 4 September 2019

# retrieve mapping statistics for all samples
# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq0" >> tmp_headers
echo "Uniq1" >> tmp_headers
echo "Uniq2" >> tmp_headers
echo "Uniq3" >> tmp_headers
echo "Uniq4" >> tmp_headers
echo "Uniq5" >> tmp_headers
echo "Uniq10" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on D2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/D2/*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsD2.tsv
rm tmp_*

# retrieve statistics for alignment on B6 mm10 primary assembly
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*exactuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq0
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*1mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq1
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*2mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*3mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq3
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*4mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq4
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*5mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq5
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenome/B6mm10_primaryassembly/*10mismatchuniqueLog.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq10

paste tmp_sample tmp_uniq0 tmp_uniq1 tmp_uniq2 tmp_uniq3 tmp_uniq4 tmp_uniq5 tmp_uniq10 > tmp_body
cat headers tmp_body > data/transcriptome/2_Mapping_STAR/OnGenome/MappingStatisticsB6mm10_primaryassembly.tsv
rm tmp_*

rm headers


# Counting reads per gene on gene subset
# 17 September 2019

# input directory
indir=data/transcriptome/2_mapping_STAR/OnGenome
# output directory
outdir=data/transcriptome/3_counting_htseq/OnGenome

# prepare subset GTF files with genes names in common between B6 and D2 transcriptome
grep -f data/transcriptome/3_counting_htseq/OnGenome/commongenes.txt references/Mus_musculus_dba2j.DBA_2J_v1.94.gtf > references/Mus_musculus_dba2j.DBA_2J_v1.94subset.gtf
grep -f data/transcriptome/3_counting_htseq/OnGenome/commongenes.txt references/Mus_musculus.GRCm38.94.gtf | grep -v havana > references/Mus_musculus.GRCm38.94subset.gtf

# generate list of commands to count reads per gene for 3 different references
for ref in B6mm10_primaryassembly D2
do
##echo $ref
	for bamfile in $indir/$ref/*_exactuniqueAligned.sortedByCoord.out.bam
	do
	sample=$(basename -s _exactuniqueAligned.sortedByCoord.out.bam $bamfile)
	##echo $sample
	if [ "$ref" == "D2" ]; then
	reffile=references/Mus_musculus_dba2j.DBA_2J_v1.94subset.gtf
	elif [ "$ref" == "B6mm10_primaryassembly" ]; then
	reffile=references/Mus_musculus.GRCm38.94subset.gtf
	fi
	echo "/home/ngobet/software/samtools-1.9/bin/samtools view -h $bamfile | htseq-count -q -s reverse -t exon -m union - $reffile 2>$outdir/$ref/$sample\_subset_htseqcount.err 1>$outdir/$ref/$sample\_subset.count" >> /mnt/nas/BXD/$outdir/CountingSubset_listcommands.sh
	done
done

# launch the counting commands in parallel
nohup parallel -j 4 < $outdir/CountingSubset_listcommands.sh > $outdir/parallel.out 2> $outdir/parallel.err &


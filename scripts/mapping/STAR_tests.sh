##############
# STAR tests #
##############

# 6 April 2020

# GOAL: test STAR adding annotation on the fly or into index

# step 1: generate genome without (a) or with (b) annotation
# step 2: map on genome with annotation on the fly (a) or on the index with annotation (b)

# working directory
cd /mnt/md0/BXD/TMP

mkdir ref map
# define variables
sample=LDB1nsd
line=DBA_2J
parental=maternal
unphased=nonrandomized
variants=indelsSNVs

# setup the genome files
cp /mnt/nas/BXD/data/PersonalizedReferencesFromGenotypes/$line\_$unphased\_$variants/*maternal.fa ref
cp /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa ref
cp /mnt/nas/BXD/data/PersonalizedReferencesFromGenotypes/$line\_$unphased\_$variants/*maternal.gtf ref

mkdir ref/star_without ref/star_with
# 1a) make genome reference STAR index without annotation
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir ref/star_without --genomeFastaFiles ref/*.fa > ref/star_without/star_without_indexing.out 2> ref/star_without/star_without_indexing.err&
# 1b) make genome reference STAR index with annotation
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir ref/star_with --genomeFastaFiles ref/*.fa --sjdbGTFfile ref/maternal.gtf > ref/star_with/star_with_indexing.out 2> ref/star_with/star_with_indexing.err&

# 2a) Run STAR alignment with annotation on the fly
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir ref/star_without --sjdbGTFfile ref/maternal.gtf --outFileNamePrefix map/onthefly_$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --quantMode GeneCounts > map/onthefly_$sample\_align_default.out 2> map/onthefly_$sample\_align_default.err&

# 2b) Run STAR alignment on the index with annotation
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir ref/star_with --outFileNamePrefix map/with_$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --quantMode GeneCounts > map/with_$sample\_align_default.out 2> map/with_$sample\_align_default.err&

# test if same output
diff -s map/onthefly_LDB1nsd_Log.final.out map/with_LDB1nsd_Log.final.out
#Yes, just the mapping times are changing.
# test time (substracting first and last time on .out files
#1a = 21 min 01 sec
#1b = 23 min 20 sec
#2a = 13 min 12 sec
#2b = 10 min 39 sec
#Total a = 34 min 13 sec
#Total b = 33 min 59 sec

# CONCLUSION: For one run, slightly better to use annotation on the fly rather than into index, but for more than one run win ~ 2 minutes by run.

#################################################################################################################################################

# 6 April 2020

# GOAL: test if HTseq count and STAR counting are really equivalent


# count with HTseq
/usr/bin/htseq-count -f bam -q -s reverse -t exon -m union map/onthefly_LDB1nsd_Aligned.sortedByCoord.out.bam ref/maternal.gtf  2>map/$sample\_htseqcount.err 1>map/$sample\_htseqcount.out&

# compare with STAR
cut -f 1,4 map/onthefly_LDB1nsd_ReadsPerGene.out.tab | sort > tmp
diff -s tmp map/LDB1nsd_htseqcount.out
# CONCLUSION: Differs only for ambigous and no features rows. However, much faster to have it integrated to STAR (because of multi-threading I guess) than HTseq alone (26 min 1 sec).
rm tmp

#################################################################################################################################################

# 6 April 2020

# GOAL: try to count genes with STAR on alignment on genome only

/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir ref/star_without --outFileNamePrefix map/genome_$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --quantMode GeneCounts > map/genome_$sample\_align_default.out 2> map/genome_$sample\_align_default.err&
#Exits with error reporting lack of gtf file.

# map on the genome with STAR and count with htseq-count if much different from STAR geneCounts with annotation
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir ref/star_without --outFileNamePrefix map/genome_$sample\_ --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 > map/genome_$sample\_align_default.out 2> map/genome_$sample\_align_default.err&
/usr/bin/htseq-count -f bam -q -s reverse -t exon -m union map/genome_LDB1nsd_Aligned.sortedByCoord.out.bam ref/maternal.gtf  2>map/$sample\_genome_htseqcount.err 1>map/$sample\_genome_htseqcount.out&
cut -f 1,4 map/onthefly_LDB1nsd_ReadsPerGene.out.tab | sort > tmp
diff -s tmp map/LDB1nsd_genome_htseqcount.out
# CONCLUSION: The numbers change a bit. Is it important?
rm tmp

/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --help
runMode

#################################################################################################################################################

# 27 April 2020

# GOAL: test if precising the standness of library makes a difference (default: Unstranded, Forward (1st read strand same as RNA), Reverse)

refdir=/mnt/nas/BXD/references
outdir=/mnt/md0/BXD/D2_withannotation/
sample=05nsd
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_D2_withannotation/ --outFileNamePrefix $outdir/$sample\testReverse_ --runThreadN 10 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --soloStrand Reverse --quantMode GeneCounts > $outdir/$sample\testReverse_align.out 2> $outdir/$sample\testReverse_align.err&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_D2_withannotation/ --outFileNamePrefix $outdir/$sample\testForward_ --runThreadN 10 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --soloStrand Forward --quantMode GeneCounts > $outdir/$sample\testForward_align.out 2> $outdir/$sample\testForward_align.err&

# CONCLUSION: It does not make any change. readStrand parameter was removed at STAR version 2.7.0e, and soloStrand can be used for single-cell samples (https://groups.google.com/forum/#!topic/rna-star/PS9p-rkGh2s)

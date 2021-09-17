# GOAL: re-run mappings and keep bam files.

# DATE: 24 November 2020 - 28 April 2021

# AUTHOR: Nastassia Gobet

#####################################################

# Mapping with reference B6 mm10 primary assembly

# variables definition
refdir=/mnt/nas/BXD/references
##samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDnsd]\{4,8\}*" | grep -e "DB[12]" -e "B6" -e "45" -e "101" -e "48" -e "100"))
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDnsd]\{4,8\}*" | grep -e "48" -e "100"))
outdir=/mnt/md0/BXD/B6mm10primaryassembly_withannotation/

mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 4 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err
	# remove unnecessary files
	rm $outdir/$sample\_SJ.out.tab
done

#####################################################

# BXD specific
refdir=/home/ngobet/projects/BXD/references
outdir=/mnt/md0/BXD/MappingEvaluation
parental=paternal
unphased=nonrandomized
variants=genotypesandimputed

mkdir /mnt/md0/BXD/MappingEvaluation/
mkdir /mnt/md0/BXD/MappingEvaluation/genotypesandimputed_withannotation_Local_0_10/
mkdir /mnt/md0/BXD/MappingEvaluation/genotypesandimputed_withoutannotation_EndToEnd_1_0/

# BXD45
cp -r data/PersonalizedReferences/BXD45_nonrandomized\_genotypesandimputed $refdir
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/45nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45nsd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/45sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/45sd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L45nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45nsd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L45sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L45sd\_SJ.out.tab&
# parental samples
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B61nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/B61nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B61nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B61nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B61nsd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB61nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/LB61nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB61nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB61nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB61nsd\_SJ.out.tab&

/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B62nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/B62nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B62nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B62nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/B62nsd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB62nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/LB62nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB62nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB62nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LB62nsd\_SJ.out.tab&

/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB1nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/DB1nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB1nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB1nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB1nsd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB1nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/LDB1nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB1nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB1nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB1nsd\_SJ.out.tab&

/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB2nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/DB2nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB2nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB2nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/DB2nsd\_SJ.out.tab&
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD45\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB2nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 16 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/LDB2nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB2nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB2nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/LDB2nsd\_SJ.out.tab&
rm -r $refdir/BXD45_nonrandomized\_genotypesandimputed

# BXD48
cp -r data/PersonalizedReferences/BXD48_nonrandomized\_genotypesandimputed $refdir
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD48\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/48nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48nsd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD48\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/48sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/48sd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD48\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L48nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48nsd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD48\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L48sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L48sd\_SJ.out.tab
rm -r $refdir/BXD48_nonrandomized\_genotypesandimputed

# BXD100
cp -r data/PersonalizedReferences/BXD100_nonrandomized\_genotypesandimputed $refdir
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD100\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/100nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100nsd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD100\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/100sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/100sd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD100\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L100nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100nsd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD100\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L100sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L100sd\_SJ.out.tab
rm -r $refdir/PersonalizedReferences/BXD100_nonrandomized\_genotypesandimputed

# BXD101
cp -r data/PersonalizedReferences/BXD101_nonrandomized\_genotypesandimputed $refdir
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD101\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/101nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101nsd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD101\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/101sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/101sd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD101\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101nsd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L101nsd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101nsd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101nsd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101nsd\_SJ.out.tab
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir /home/ngobet/projects/BXD/references/BXD101\_nonrandomized\_genotypesandimputed/star_paternal\_withannotation --outFileNamePrefix /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101sd\_ --outSAMtype BAM SortedByCoordinate --runThreadN 1 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L101sd.fastq.gz --quantMode GeneCounts --scoreDelOpen -2 --scoreInsOpen -2 --alignIntronMax 0 --alignEndsType Local --outFilterMismatchNmax 10 > /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101sd\_align_default.out 2> /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101sd\_align.err; rm /mnt/md0/BXD/MappingEvaluation/genotypesandimputed\_withannotation\_Local\_0\_10/L101sd\_SJ.out.tab
rm -r $refdir/PersonalizedReferences/BXD101_nonrandomized\_genotypesandimputed

#####################################################

# Mapping with reference D2 (mm10 modified with D2-specific indels and SNPs from dbSNP), no modification for heterozygous sites

# variables definition
refdir=/mnt/nas/BXD/data/PersonalizedReferences
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDnsd]\{4,8\}*" | grep -e "DB[12]" -e "B6"))
outdir=/mnt/md0/BXD/DBA_2J_nonrandomized_indelsSNVs_paternal_withannotation/

mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/DBA_2J_nonrandomized_indelsSNVs/star_paternal_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 4 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err
done

# mapping BXD81 samples on BXD81-specific and GRCm38

# variables definition
refdir=/mnt/nas/BXD/data/PersonalizedReferences
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L]\{0,1\}81nsd"))
outdir=/mnt/md0/BXD/BXD81_nonrandomized_genotypesandimputed

mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/BXD81_nonrandomized_genotypesandimputed/star_paternal_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 4 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err&
done

# mapping BXD81 samples on GRCm38

# variables definition
outdir=/mnt/md0/BXD/B6mm10primaryassembly_withannotation/
refdir=/mnt/nas/BXD/references

mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 4 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err&
done

cp /mnt/md0/BXD/B6mm10primaryassembly_withannotation/*81* tmpBAM/B6/&
mkdir tmpBAM/BXD81
cp /mnt/md0/BXD/BXD81_nonrandomized_genotypesandimputed/* tmpBAM/BXD81/&

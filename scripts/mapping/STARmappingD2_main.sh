# STAR mapping on D2 assembly

# GOAL: re-run mapping on D2 assembly with default parameters and 1-pass.

# DATE: 24 April 2020 - 10 March 2021

# build reference index for STAR with annotation
refdir=/mnt/nas/BXD/references
mkdir $refdir/genome/star2.7.0e_D2_withannotation/
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/genome/star2.7.0e_D2_withannotation/ --genomeFastaFiles $refdir/genome/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa --sjdbGTFfile $refdir/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.94.gtf --outFileNamePrefix $refdir/genome/star2.7.0e_D2_withannotation/star2.7.0e_D2_withannotation_ > $refdir/genome/star2.7.0e_D2_withannotation/indexing.out 2> $refdir/genome/star2.7.0e_D2_withannotation/indexing.err
##analysis/scripts/mapping/prepareReference.sh

# run alignment on all samples
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDtnsd]\{4,8\}*"))
outdir=/mnt/md0/BXD/D2_withannotation/
mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_D2_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err
	# remove unnecessary files
	rm $outdir/$sample\_Aligned.out.bam
	rm $outdir/$sample\_SJ.out.tab
done

# retrieve uniquely mapped percentage for alignment on genome D2 assembly with annotation

# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on D2
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenomeAndAnnotation/D2_withannotation/*Log.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenomeAndAnnotation/D2_withannotation/*Log.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq

paste tmp_sample tmp_uniq* > tmp_body
tr -s '\t' < tmp_body > tmp_body2
cat headers tmp_body2 > data/transcriptome/2_Mapping_STAR/OnGenomeAndAnnotation/MappingStatisticsD2assembly.tsv
rm tmp_* headers

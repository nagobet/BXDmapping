# STAR mapping on B6 mm9

# GOAL: re-run mapping on mm9 non-masked genome+annotation with STAR default parameters (for Figure 2A).

# DATE: 24 April 2020

# download reference sequence
##curl -O ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz
# download transcriptome annotation on reference genome
##curl -O ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz

# build reference index for STAR with annotation
refdir=/mnt/nas/BXD/references
##mkdir $refdir/genome/star2.7.0e_B6mm9_withannotation/
##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/genome/star2.7.0e_B6mm9_withannotation/ --genomeFastaFiles $refdir/genome/Mus_musculus.NCBIM37.67.dna.toplevel.fa --sjdbGTFfile $refdir/transcriptome/Mus_musculus.NCBIM37.67.gtf --outFileNamePrefix $refdir/genome/star2.7.0e_B6mm9_withannotation/star2.7.0e_B6mm9_withannotation_ > $refdir/genome/star2.7.0e_B6mm9_withannotation/indexing.out 2> $refdir/genome/star2.7.0e_B6mm9_withannotation/indexing.err
##analysis/scripts/mapping/prepareReference.sh

# run alignment on all samples
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDtnsd]\{4,8\}*"))
outdir=/mnt/md0/BXD/B6mm9_withannotation/
mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_B6mm9_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err
	# remove unnecessary files
	rm $outdir/$sample\_Aligned.out.bam
	rm $outdir/$sample\_SJ.out.tab
done

# Run alignment

# GOAL: run STAR alignment with parameters and options
# options: reference line, unphased, variants, parental
# -s (sample name)
# -o (output directory)
# -l (reference line) -a (annotation) -u (unphased behavior)  -p ('parental' sequence according to vcf2diploid) -v (variants used to customize reference)
# to implement: -r (reference index directory) that can be passed instead of -l, -a, -u, -p, -v

# retrieve short flag arguments values
while getopts s:o:r:l:a:u:p:v: option
do
case "${option}" in
	s) sample=${OPTARG};;
	o) outputdir=${OPTARG};;
	r) refdir=${OPTARG};;
	l) line=${OPTARG};;
	a) annotation=${OPTARG};;
	u) unphased=${OPTARG};;
	p) parental=${OPTARG};;
	v) variants=${OPTARG};;
esac
done
echo "The sample is $sample."
echo "The output directory is $outputdir."
echo "The reference directory is $refdir/$line\_genome$annotation\_$unphased\_$variants\_$parental."
echo "(line = $line, annotation = $annotation, unphased behavior = $unphased, parental = $parental, variants = $variants)"

counting="-"
if [ $annotation = "withannotation" ]
then
counting="GeneCounts"
fi
echo $counting

# Run STAR alignment
mkdir $outputdir/$line\_genome$annotation\_$unphased\_$variants\_$parental/
##samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L]\{0,1\}[0-9BDnsd]\{3,6\}*"))
# restrictive (exact matches) on genome
##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental --outFileNamePrefix $outputdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 > $outputdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_align_restrictive.out 2> $outputdir/$line\_genome_$unphased\_$variants\_$parental/$sample\_align_restrictive.err
# defaults on genome + annotation
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_$annotation --outFileNamePrefix $outputdir/$line\_genome$annotation\_$unphased\_$variants\_$parental/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --quantMode $counting > $outputdir/$line\_genome$annotation\_$unphased\_$variants\_$parental/$sample\_align.out 2> $outputdir/$line\_genome$annotation\_$unphased\_$variants\_$parental/$sample\_align.err
# remove unnecessary files
rm $outputdir/$line\_genome$annotation\_$unphased\_$variants\_$parental/$sample\_Aligned.out.bam
rm $outputdir/$line\_genome$annotation\_$unphased\_$variants\_$parental/$sample\_SJ.out.tab

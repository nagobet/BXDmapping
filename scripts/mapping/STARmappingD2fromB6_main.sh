# STAR mapping on D2 from B6

# GOAL: re-run mapping on D2 reference (customized from B6 mm10 primary assembly, using indels and SNVs from dbSNP)  with default parameters and 1-pass.

# DATE: 11 May 2020

# Variables definition
refdir=/mnt/nas/BXD/data/PersonalizedReferences/DBA_2J_nonrandomized_indelsSNVs
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -v -e "05" -e "29" -e "32" -e "63" | grep -o "[L0-9BDnsd]\{4,8\}*"))
outdir=/mnt/md0/BXD/DBA_2J_nonrandomized_indelsSNVs

mkdir $outdir
echo "#DATE:" `date` > /mnt/nas/BXD/analysis/scripts/mapping/STARmappingD2fromB6_listcommands.sh

# Restrictive mapping settings
# run alignment on all samples
for sample in "${samples[@]}"
do
	##echo $sample
	# alignment on genome + annotation
	echo "/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/star_paternal_withoutannotation --outFileNamePrefix $outdir/$sample\_restrictive_ --runThreadN 2 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd Log --scoreDelOpen -40 --scoreInsOpen -40 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 0 > $outdir/$sample\_restrictive_align.out 2> $outdir/$sample\_restrictive_align.err; htseq-count -f bam -q -s reverse -t exon -m union $outdir/$sample\_restrictive_Aligned.out.bam $refdir/paternal.gtf 2> $outdir/$sample\_htseqcount.err 1> $outdir/$sample\_restrictive_ReadsPerGene.out.tab; rm $outdir/$sample\_restrictive*.bam; rm $outdir/$sample\_restrictive*SJ.out.tab" >> /mnt/nas/BXD/analysis/scripts/mapping/STARmappingD2fromB6_listcommands.sh
done

# Permissive (default) mapping settings
# run alignment on all samples
for sample in "${samples[@]}"
do
	##echo $sample
	# alignment on genome + annotation
	echo "/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/star_paternal_withannotation --outFileNamePrefix $outdir/$sample\_permissive_ --runThreadN 2 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --quantMode GeneCounts > $outdir/$sample\_permissive_align.out 2> $outdir/$sample\_permissive_align.err; rm $outdir/$sample\_permissive_Aligned.out.bam;	rm $outdir/$sample\_permissive_SJ.out.tab" >> /mnt/nas/BXD/analysis/scripts/mapping/STARmappingD2fromB6_listcommands.sh
done

# launching jobs
nohup parallel -j 4 --delay 1 < /mnt/nas/BXD/analysis/scripts/mapping/STARmappingD2fromB6_listcommands.sh

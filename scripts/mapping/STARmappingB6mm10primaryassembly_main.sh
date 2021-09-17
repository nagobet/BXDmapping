# STAR mapping on B6 mm10 

# GOAL: re-run mapping on mm10 primary assembly

# DATE: 24 April 2020 - 10 March 2021

# build reference index for STAR with annotation
refdir=/mnt/nas/BXD/references
mkdir $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/
/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/ --genomeFastaFiles $refdir/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa --sjdbGTFfile $refdir/transcriptome/Mus_musculus.GRCm38.94.gtf --outFileNamePrefix $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/star2.7.0e_B6mm10primaryassembly_withannotation_ > $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/indexing.out 2> $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/indexing.err

# GTFtools only consider genes on conventional genes (help from: https://genomespot.blogspot.com/2019/01/using-gtf-tools-to-get-gene-lengths.html)
grep -v '#' references/transcriptome/Mus_musculus.GRCm38.94.gtf | awk '{OFS="\t"} $1=1' > tmp.gtf
/home/ngobet/software/GTFtools_0.6.5/gtftools.py -l references/transcriptome/genelength.txt tmp.gtf
rm tmp.gtf

# run alignment on all samples
samples=($(ls /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*.fastq.gz | grep -o "[L0-9BDtnsd]\{4,8\}*"))
outdir=/mnt/md0/BXD/B6mm10primaryassembly_withannotation/
mkdir $outdir
for sample in "${samples[@]}"
do
	echo $sample
	# alignment on genome + annotation
	/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --genomeDir $refdir/genome/star2.7.0e_B6mm10primaryassembly_withannotation/ --outFileNamePrefix $outdir/$sample\_ --runThreadN 18 --readFilesCommand zcat --readFilesIn /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/$sample.fastq.gz --outSAMtype BAM Unsorted --outStd Log --quantMode GeneCounts > $outdir/$sample\_align.out 2> $outdir/$sample\_align.err
	# remove unnecessary files
	rm $outdir/$sample\_Aligned.out.bam
	rm $outdir/$sample\_SJ.out.tab
done

# group counts
analysis/scripts/mapping/MergeCount.sh data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/ _ReadsPerGene.out.tab
mv data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/Summary_ReadsPerGene.out.tab data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/SummaryAll_ReadsPerGene.out.tab
# filter out parental and F1 samples
module load R/3.4.2
Rscript analysis/scripts/mapping/filterSamplesCounts.R

# normalize counts CPM
Rscript analysis/scripts/mapping/normalizeGeneCountsCPM.R data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/
# normalize counts TPM
Rscript analysis/scripts/mapping/normalizeGeneCountsTPM_B6mm10.R data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/

# group and filter in restrictive mapping setting
# group counts
analysis/scripts/mapping/MergeCount.sh data/transcriptome/2_mapping_STAR/OnGenome/B6mm10_primaryassembly/ _ReadsPerGene.out.tab
mv data/transcriptome/2_mapping_STAR/OnGenome/B6mm10_primaryassembly/Summary_ReadsPerGene.out.tab data/transcriptome/2_mapping_STAR/OnGenome/B6mm10_primaryassembly/SummaryAll_ReadsPerGene.out.tab
# filter out parental and F1 samples
module load R/3.4.2
Rscript analysis/scripts/mapping/filterSamplesCounts_B6mm10_restrictive.R


########################################################################################
# To run from Wally cluster

# copy expression files on cluster from 3009
codedir=/users/ngobet/scripts
tmpdir=bxdmap_scratchdirectory2/data/eQTL_B6mm10primaryassembly_withannotation
mkdir $tmpdir
scp -r ngobet@pccig3009.unil.ch:/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/TMMnormalized_log2*PM_*.tab $tmpdir

# eQTL detection
setting=$(ls -d $tmpdir | cut -d "/" -f 3)
for tissue in Cortex Liver
do
	echo $tissue
	for condition in NSD SD
	do
		echo $condition
		# transform phenotypes (gene expression) to bed format (UCSC format)
		scripts/transformGenePhenotypesToBedUCSC.py $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.tab $refdir/genotypes/BXDGenotypes.geno $refdir/transcriptome_annotation/GenePosition.txt $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.bed > $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition\_transformGenePhenotypesToBedUCSC.out 2> $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition\_transformGenePhenotypesToBedUCSC.err
		scripts/transformGenePhenotypesToBedUCSC.py $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.tab $refdir/genotypes/BXDGenotypes.geno $refdir/transcriptome_annotation/GenePosition.txt $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.bed > $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition\_transformGenePhenotypesToBedUCSC.out 2> $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition\_transformGenePhenotypesToBedUCSC.err
		# compress and index phenotypes
		bgzip -f $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.bed && tabix -p bed $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.bed.gz
		bgzip -f $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.bed && tabix -p bed $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.bed.gz

		# prepare commands
		fastQTL --vcf $refdir/genotypes/BXDGenotypes.vcf.gz --bed $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.bed.gz --out $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition --commands 25 $codedir/eQTL_$setting\_$tissue\_$condition\_CPM_listcommands.sh --window 2e6 --permute 1000 --seed 1
		fastQTL --vcf $refdir/genotypes/BXDGenotypes.vcf.gz --bed $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.bed.gz --out $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition --commands 25 $codedir/eQTL_$setting\_$tissue\_$condition\_TPM_listcommands.sh --window 2e6 --permute 1000 --seed 1
		# run commands
		bash $codedir/eQTL_$setting\_$tissue\_$condition\_CPM_listcommands.sh > $tmpdir/eQTL_$setting\_$tissue\_$condition\_CPM_listcommands.out 2> $tmpdir/eQTL_$setting\_$tissue\_$condition\_CPM_listcommands.err
		bash $codedir/eQTL_$setting\_$tissue\_$condition\_TPM_listcommands.sh > $tmpdir/eQTL_$setting\_$tissue\_$condition\_TPM_listcommands.out 2> $tmpdir/eQTL_$setting\_$tissue\_$condition\_TPM_listcommands.err
		
		# post-processing
		# group results from different regions
		cat $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.chr* > $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.txt
		cat $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.chr* > $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.txt
		# clean up unneeded files
		rm $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition.chr* $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition.chr* $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition
		# calculate q-values
		Rscript $codedir/correctMultiPhenotypeseQTL.R $tmpdir/TMMnormalized_log2CPM\_$tissue\_$condition
		Rscript $codedir/correctMultiPhenotypeseQTL.R $tmpdir/TMMnormalized_log2TPM\_$tissue\_$condition
	done
done

# remove expression files on cluster
rm $tmpdir/TMMnormalized_log2*PM_*.bed* $tmpdir/TMMnormalized_log2*PM_*.tab

# copy expression files on 3009
scp -r $tmpdir ngobet@pccig3009.unil.ch:/mnt/nas/BXD/data/transcriptome/$setting

# remove eQTL files on cluster
rm -r $tmpdir

########################################################################################

# retrieve uniquely mapped percentage for alignment on genome B6 mm10 with annotation

# prepare headers
echo "Sample_Name" > tmp_headers
echo "Uniq" >> tmp_headers
paste -s tmp_headers > headers

# retrieve statistics for alignment on B6 mm10
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/*Log.final.out | cut -d / -f 6 | cut -d _ -f 1 > tmp_sample
grep "Uniquely mapped reads %" data/transcriptome/2_Mapping_STAR/OnGenomeAndAnnotation/B6mm10primaryassembly_withannotation/*Log.final.out | cut -d / -f 6 | cut -d "|" -f 2 | cut -d "%" -f 1 > tmp_uniq

paste tmp_sample tmp_uniq* > tmp_body
tr -s '\t' < tmp_body > tmp_body2
cat headers tmp_body2 > data/transcriptome/2_Mapping_STAR/OnGenomeAndAnnotation/MappingStatisticsB6mm10primaryassembly.tsv
rm tmp_* headers

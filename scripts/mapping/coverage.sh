
# generate commands to have the coverage
##for ref in D2 B6mm10 D2_majorchromosomes B6mm10_majorchromosomes; do
for ref in D2 B6mm10_primaryassembly; do
	##echo $ref
	# take only the parental nsd samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $sample
	echo "/home/ngobet/software/bedtools2/bin/genomeCoverageBed -dz -ibam data/transcriptome/2_mapping_STAR/OnGenome/$ref/$sample\_exactuniqueAligned.sortedByCoord.out.bam > data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_exactunique.cov" >> scripts/listcommands_coverage.sh
	done
done

nohup parallel -j 16 < scripts/listcommands_coverage.sh 1> cov_parallel.out 2> cov_parallel.err &
#One job is taking execessively long to run: DB1nsd on B6mm10_primaryassembly (I don't know why)

# split by chromosomes (to have smaller files)
##for ref in D2 B6mm10 D2_majorchromosomes B6mm10_majorchromosomes; do
##for ref in D2_majorchromosomes; do
for ref in D2 B6mm10_primaryassembly; do
	for myfile in data/transcriptome/3_coverage_bedtools/OnGenome/$ref/*[BD]*[12]nsd_exactunique.cov; do
	sample=$(basename -s _exactunique.cov $myfile)
	##echo $myfile
	##echo $sample
	awk '{print >> "data/transcriptome/3_coverage_bedtools/OnGenome/'$ref'/""'$sample'_"$1".cov"; close("data/transcriptome/3_coverage_bedtools/OnGenome/'$ref'/""'$sample'_"$1".cov")}' $myfile &
	done
done
# works but long. ~1KB/s so 31 days for finishing. :(
#help from: https://stackoverflow.com/questions/16635396/split-large-file-according-to-value-in-single-column-awk


########### test ###########
# take the coverage directly by chrosomomes
##split -l 1 -d -a 2 references/genome/star2.7.0e_D2_majorchromosomes/chrNameLength.txt chrNameLength
##split -l 1 -d -a 2 references/genome/star2.7.0e_B6mm10_primaryassembly/chrNameLength.txt chrNameLength

##awk '{print >> "chr"$1".txt"}' references/genome/star2.7.0e_D2_majorchromosomes/chrNameLength.txt
head -n 20 references/genome/star2.7.0e_D2/chrNameLength.txt | awk '{print >> "chr"$1".txt"}'
##head -n 22 references/genome/star2.7.0e_B6mm10_primaryassembly/chrNameLength.txt | awk '{print >> "chr"$1".txt"}'

ref="D2"
ref="B6mm10_primaryassembly"
##echo $ref
for mychr in *.txt; do
chr=$(basename -s .txt $mychr)
	# take only the parental nsd samples
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz; do
	sample=$(basename -s .fastq.gz $fqfile)
	##echo $sample
	echo "/home/ngobet/software/bedtools2/bin/genomeCoverageBed -dz -ibam data/transcriptome/2_mapping_STAR/OnGenome/$ref/$sample\_exactuniqueAligned.sortedByCoord.out.bam -g $mychr > data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_exactunique.cov" >> scripts/listcommands_coverage.sh
	done
done
rm chr*.txt

### NOT WORKING "WARNING: Genome (-g) files are ignored when BAM input is provided."

########### test 2 ################
# find a faster way to split by chromosome

##for ref in D2 B6mm10 D2_majorchromosomes B6mm10_majorchromosomes
##for ref in B6mm10 D2_majorchromosomes B6mm10_majorchromosomes
for ref in B6mm10_primaryassembly D2
do
	echo $ref
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz
	do
	sample=$(basename -s .fastq.gz $fqfile)
	echo $sample
	for mychr in {1..19} X
	do
	echo $mychr
	grep -n -w "^$mychr" data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_exactunique.cov | head -n 1 >> data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_lines_start
	grep -n -w "^$mychr" data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_exactunique.cov | tail -n 1 >> data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_lines_end
	done
	done
done
# takes quite a while (~ 1h30)

##for ref in D2 B6mm10 D2_majorchromosomes B6mm10_majorchromosomes
##for ref in D2 B6mm10 D2_majorchromosomes
for ref in B6mm10_primaryassembly D2
do
	echo $ref
	for fqfile in /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/*[BD]*[12]nsd.fastq.gz
	do
	sample=$(basename -s .fastq.gz $fqfile)
	echo $sample
	for mychr in {1..19} 0
	do
	echo $mychr
	cut -f 1 -d : data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_lines_start > tmpstart
	cut -f 2 -d : data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_lines_start | cut -f 1 > keys
	cut -f 1 -d : data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_lines_end > tmpend
	paste keys tmpstart tmpend > tmp1
	# take X chromosome and put in on the beginning of the file
	tail -n 1 tmp1 > tmp
	head -n 19 tmp1 >> tmp
	# make 2 arrays with start and end lines for each chromosome. index is chromosome number, except for X chromosome which has index 0.
	while read c1 c2 c3 leftovers;do starts[$c1]=$c2; ends[$c1]=$c3; done < tmp
	##sed -n "${starts[$mychr]},${ends[$mychr]}p;$[${ends[$mychr]}+1]q" data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_exactunique.cov > data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_$mychr.cov
	sed -n "${starts[$mychr]},${ends[$mychr]}p;$[${ends[$mychr]}+1]q" data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_exactunique.cov > data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_$mychr.cov
	done
	rm tmp tmp1 keys tmpstart tmpend
	##mv data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_0.cov data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_X.cov
	mv data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_0.cov data/transcriptome/3_coverage_bedtools/OnGenome/$ref/$sample\_X.cov
	done
done
# ~1h

##################################################

# 29 August 2019

# check if annotation is in agreement with coverage on genome

# get genes (column V4 and V5 are positions)
grep -w "gene" Mus_musculus_dba2j.DBA_2J_v1.94.gtf > geneD2.tab
grep -w "gene" Mus_musculus.GRCm38.94.gtf > geneB6mm10.tab
# get exons
grep -w "exon" Mus_musculus_dba2j.DBA_2J_v1.94.gtf > exonD2.tab
grep -w "exon" Mus_musculus.GRCm38.94.gtf > exonB6mm10.tab

# convert sam to bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_exactuniqueAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/B61nsd_exactuniqueAligned.sortedByCoord.out.bam
/home/ngobet/software/samtools-1.9/bin/samtools view -o /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/LB61nsd_exactuniqueAligned.sortedByCoord.out.sam /mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnGenome/D2/LB61nsd_exactuniqueAligned.sortedByCoord.out.bam

# count how many alignments fall into exons
# bedtools intersect

# convert bam to bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnGenome/D2/LB61nsd_exactuniqueAligned.sortedByCoord.out.bam -cigar > LB61nsd_onD2_0.bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnGenome/D2/LB62nsd_exactuniqueAligned.sortedByCoord.out.bam -cigar > LB62nsd_onD2_0.bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnGenome/D2/LDB1nsd_exactuniqueAligned.sortedByCoord.out.bam -cigar > LDB1nsd_onD2_0.bed &

# overlap with genes
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LB61nsd_onD2_0.bed -b references/geneD2.tab > overlap_genes.bed &
cut -f 4 overlap_genes.bed  > reads_overlap_genes.bed 
uniq reads_overlap_genes.bed  | wc -l
#11299959
wc -l LB61nsd_onD2_0.bed
#12422079
# So 91% of reads mapped overlap with genes

# overlap with exons
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LB61nsd_onD2_0.bed -b references/exonD2.tab > overlap_exons.bed &
cut -f 4 overlap_exons.bed > reads_overlap_exons.bed
uniq reads_overlap_exons.bed | wc -l
#10835429
# So 87% of reads mapped overlap with exons

# replicate
# overlap with genes
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LB62nsd_onD2_0.bed -b references/geneD2.tab > overlap_genes.bed &
cut -f 4 overlap_genes.bed | uniq > reads_overlap_genes.bed 
wc -l reads_overlap_genes.bed
#10667461
wc -l LB62nsd_onD2_0.bed
#11706280
# So 91% of reads mapped overlap with genes

# overlap with exons
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LB62nsd_onD2_0.bed -b references/exonD2.tab > overlap_exons.bed &
cut -f 4 overlap_exons.bed | uniq > reads_overlap_exons.bed
wc -l reads_overlap_exons.bed
#10215649
# So 87% of reads mapped overlap with exons

# on DB sample
# overlap with genes
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LDB1nsd_onD2_0.bed -b references/geneD2.tab > overlap_genes.bed &
cut -f 4 overlap_genes.bed > reads_overlap_genes.bed &
uniq reads_overlap_genes.bed | wc -l
#12164107
wc -l LDB1nsd_onD2_0.bed
#13246698
# So 92% of reads mapped overlap with genes

# overlap with exons
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LDB1nsd_onD2_0.bed -b references/exonD2.tab > overlap_exons.bed &
cut -f 4 overlap_exons.bed | uniq > reads_overlap_exons.bed &
wc -l reads_overlap_exons.bed
#11562805 
# So 87% of reads mapped overlap with exons


# check on transcriptome
# run alignment on D2 transcriptome
module add Alignment/STAR/2.7.0e
STAR --genomeDir references/star2.7.0e_transcriptome_D2_old --outFilterIntronMotifs RemoveNoncanonicalUnannotated --twopassMode Basic --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/LB61nsd --runThreadN 8 --readFilesCommand zcat --readFilesIn data/transcriptome/RNAseq/filtered/LB61nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMultimapNmax 1 &
STAR --genomeDir references/star2.7.0e_transcriptome_D2_old --outFilterIntronMotifs RemoveNoncanonicalUnannotated --twopassMode Basic --outFileNamePrefix /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/LDB1nsd --runThreadN 8 --readFilesCommand zcat --readFilesIn data/transcriptome/RNAseq/filtered/LDB1nsd.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 50000000000 --outFilterMultimapNmax 1 &
mv /home/ngobet/projects/BXD/data/TMP/OnTranscriptome/D2/L* data/transcriptome/2_mapping_STAR/OnTranscriptome/D2/ &

# convert bam to bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnTranscriptome/D2/LB61nsdAligned.sortedByCoord.out.bam -cigar > LB61nsd_onD2_0_transcriptome.bed
/home/ngobet/software/bedtools2/bin/bamToBed -i data/transcriptome/2_mapping_STAR/OnTranscriptome/D2/LDB1nsdAligned.sortedByCoord.out.bam -cigar > LDB1nsd_onD2_0_transcriptome.bed

# LB61nsd
# overlap with genes
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LB61nsd_onD2_0_transcriptome.bed -b references/geneD2.tab > overlap_genes.bed &
cut -f 4 overlap_genes.bed | uniq > reads_overlap_genes.bed 
wc -l reads_overlap_genes.bed
#30669154
wc -l LB61nsd_onD2_0_transcriptome.bed
#35482044
# So 86% of reads mapped overlap with genes

# overlap with exons
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LB61nsd_onD2_0_transcriptome.bed -b references/exonD2.tab > overlap_exons.bed &
cut -f 4 overlap_exons.bed | uniq > reads_overlap_exons.bed
wc -l reads_overlap_exons.bed
#29196855
# So 82% of reads mapped overlap with exons

# LDB1nsd
# overlap with genes
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LDB1nsd_onD2_0_transcriptome.bed -b references/geneD2.tab > overlap_genes.bed &
cut -f 4 overlap_genes.bed > reads_overlap_genes.bed 
uniq reads_overlap_genes.bed | wc -l
#28101951
wc -l LDB1nsd_onD2_0_transcriptome.bed
#32591940
# So 86% of reads mapped overlap with genes

# overlap with exons
/home/ngobet/software/bedtools2/bin/intersectBed -wo -a LDB1nsd_onD2_0_transcriptome.bed -b references/exonD2.tab > overlap_exons.bed &
cut -f 4 overlap_exons.bed > reads_overlap_exons.bed
uniq reads_overlap_exons.bed | wc -l
#26766608
# So 82% of reads mapped overlap with exons


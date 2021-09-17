# Prepare Reference

# GOAL: prepare reference for STAR mapping

# define frequently used directories
refdir=/mnt/nas/BXD/data/PersonalizedReferences
aligndir=/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/OnPersonalizedRef
tmpdir=/mnt/md0/BXD/TMP

##lines=($(cut -f 6 data/ConvertLineNames.tsv | grep "BXD" | awk 'NR > 4 { print }'))
##for line in "${lines[@]}"
for line in B6D2F1
do
	echo $line
	##linenumber=$(grep -o [0-9]* <<< $line)
	##echo $linenumber
	# build personalized genomes
	##export JAVA_HOME=/home/ngobet/software/java/jdk-13.0.2
	##export PATH=$JAVA_HOME/bin:$PATH
	# unphased variant randomized or not
	##mkdir $refdir/$line\_randomized_genotypes $refdir/$line\_nonrandomized_genotypes
	##cd $refdir/$line\_randomized_genotypes
	##java -jar /home/ngobet/software/vcf2diploid-master/vcf2diploid.jar \
	## -id $line \
	## -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	## -vcf /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf > vcf2diploid.out 2> vcf2diploid.err
	## cd $refdir/$line\_nonrandomized_genotypes
	## java -jar /home/ngobet/software/vcf2diploid-masterNotRandomized/vcf2diploid.jar \
	## -id $line \
	## -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	## -vcf /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf > vcf2diploid.out 2> vcf2diploid.err
	##for parental in maternal paternal
	for parental in paternal
	do
	##for unphased in nonrandomized randomized
	for unphased in nonrandomized
	do
	echo $parental\_$unphased\_genotypes
	##mkdir $refdir/$line\_$unphased\_genotypes/star_$parental $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation
	# liftover transcriptome annotation
	##/home/ngobet/software/liftOver -gff /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf $refdir/$line\_$unphased\_genotypes/$parental\.chain $refdir/$line\_$unphased\_genotypes/$parental\.gtf $refdir/$line\_$unphased\_genotypes/$parental\_unlifted.gtf
	# build references index for STAR
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_genotypes/star_$parental --genomeFastaFiles $refdir/$line\_$unphased\_genotypes/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --outFileNamePrefix $refdir/$line\_$unphased\_genotypes/star_$parental/star_$parental\_ > $refdir/$line\_$unphased\_genotypes/star_$parental/indexing.out 2> $refdir/$line\_$unphased\_genotypes/star_$parental/indexing.err
	##/software/Alignment/STAR/STAR-2.7.0e/bin/STAR --runThreadN 18 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation --genomeFastaFiles $refdir/$line\_$unphased\_genotypes/*$parental\.fa /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --sjdbGTFfile $refdir/$line\_$unphased\_genotypes/$parental.gtf --outFileNamePrefix $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation/star_$parental\_withannotation_ > $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation/indexing.out 2> $refdir/$line\_$unphased\_genotypes/star_$parental\_withannotation/indexing.err
	
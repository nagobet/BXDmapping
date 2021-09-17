#!/bin/bash

###############################################################################

# CONTENT header

# STAR mapping optimization

# GOAL: run mapping various STAR mapping parameters to identify optimal parameters

# DATE: 25 May 2020 -26 June 2020

# INFO:
# Steps from reads of samples (.fastq.gz) to cis eQTL detection
# Some steps are run locally on 3009 computer (see MappingOptimization_main.sh)

###############################################################################

# Part 1: Mapping

# directories definitions
refdir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/references
indir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/reads
mapdir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/mapping
eQTLdir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/eQTL
codedir=/users/ngobet/scripts
# softwares definitions
STARpath=/users/ngobet/softwares/STAR-2.7.0e/bin/Linux_x86_64
liftOverpath=/users/ngobet/softwares
vcf2diploidpath=/users/ngobet/softwares/vcf2diploid-masterNotRandomized
# variables definitions
parental=paternal
unphased=nonrandomized
variants=genotypesandimputed
deletion=2
insertion=2

# Prepare mapping commands
for line in $(cat data/lines.txt)
do
	echo "#line: $line"
	echo "# commands" > $codedir/MappingOptimization_$line\_listcommands.sh
	# Run STAR alignment
	for annotation in withoutannotation withannotation
	do
		echo "##$annotation"
		for trimming in Local EndToEnd
		do
			echo "###trimming: $trimming"
			for intronMax in 1 0
			do
				echo "####intronMax: $intronMax"
				for mismatches in 0 1 2 3 10
				do
					echo "#####mismatches: $mismatches"
					mkdir -p $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/
					chmod 777 $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/
					linenumber=$(grep -o [0-9]* <<< $line)
					samples=($(ls $indir/*.fastq.gz | grep -o "[L0-9D]\{2,4\}[n]\{0,1\}sd*" | grep "$linenumber" | sort | uniq))
					for sample in "${samples[@]}"
					do
						# mapping command (alignment with STAR + gene counting)
						# mapping without transcriptome annotation
						if [ $annotation = "withoutannotation" ]
						then
							echo "$STARpath/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_$annotation --outFileNamePrefix $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_ --outSAMmode Full --outStd SAM --runThreadN 2 --readFilesCommand zcat --readFilesIn $indir/$sample.fastq.gz --scoreDelOpen -$deletion --scoreInsOpen -$insertion --alignIntronMax $intronMax --alignEndsType $trimming --outFilterMismatchNmax $mismatches 2> $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_align.err | htseq-count -f bam -q -s reverse -t exon -m union - $refdir/$line\_$unphased\_$variants/$parental.gtf 2> $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_htseqcount.err 1> $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_ReadsPerGene.out.tab; rm $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_SJ.out.tab" >> $codedir/MappingOptimization_$line\_listcommands.sh
						# mapping with transcriptome annotation
						elif [ $annotation = "withannotation" ]
						then
							echo "$STARpath/STAR --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_$annotation --outFileNamePrefix $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_ --outSAMtype BAM Unsorted --runThreadN 2 --readFilesCommand zcat --readFilesIn $indir/$sample.fastq.gz --quantMode GeneCounts --scoreDelOpen -$deletion --scoreInsOpen -$insertion --alignIntronMax $intronMax --alignEndsType $trimming --outFilterMismatchNmax $mismatches > $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_align_default.out 2> $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_align.err; rm $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_*.bam; rm $mapdir/$variants\_$annotation\_$trimming\_$intronMax\_$mismatches/$sample\_SJ.out.tab" >> $codedir/MappingOptimization_$line\_listcommands.sh
						fi
					done #end samples loop
				done #end mismatches loop
			done #end intronMax loop
		done #end trimming loop
	done #end annotation loop
done #end lines loop

# Run mapping commands by line (reference)
##for line in $(cat data/lines.txt)
for line in $(cat data/lines.txt)
do
	echo "#$line"
	# # Impute genotypes for BXD line
	# # retrieve header
	# grep "#CHROM" data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# # adapt header to BXD line
	# sed -i "s/DBA_2J/$line/g" data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# # impute genotypes for BXD line
	# /home/ngobet/software/bedtools2/bin/intersectBed -a data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf -b data/genome/$line\_D2blocks.bed -wa >> data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf
	# /home/ngobet/software/bedtools2/bin/intersectBed -a data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf -b data/genome/$line\_D2blocks.bed -wa >> data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf

	# # Build personalized genomes
	# mkdir $refdir/$line\_$unphased\_$variants
	# cd $refdir/$line\_$unphased\_$variants
	# java -jar $vcf2diploidpath/vcf2diploid.jar \
	 # -id $line \
	 # -chr /mnt/nas/BXD/references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	 # -vcf /mnt/nas/BXD/data/genome/ImputedGenotypes/$line\_imputedgenotypes.vcf > $refdir/$line\_$unphased\_$variants/vcf2diploid.out 2> $refdir/$line\_$unphased\_$variants/vcf2diploid.err
	# cd /mnt/nas/BXD/
	# mkdir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation
	# # liftover transcriptome annotation
	# $liftOverpath/liftOver -gff $refdir/transcriptome_annotation/Mus_musculus.GRCm38.94.gtf $refdir/$line\_$unphased\_$variants/$parental\.chain $refdir/$line\_$unphased\_$variants/$parental\.gtf $refdir/$line\_$unphased\_$variants/$parental\_unlifted.gtf
	# # build references index for STAR
	# $STARpath/STAR --runThreadN 8 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa $refdir/genome_sequences/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/star_$parental\_withoutannotation_\_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withoutannotation/indexing.err
	# $STARpath/STAR --runThreadN 8 --limitGenomeGenerateRAM 90000000000 --runMode genomeGenerate --genomeDir $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation --genomeFastaFiles $refdir/$line\_$unphased\_$variants/*$parental\.fa $refdir/genome_sequences/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa --sjdbGTFfile $refdir/$line\_$unphased\_$variants/$parental.gtf --outFileNamePrefix $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/star_$parental\_withannotation_ > $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.out 2> $refdir/$line\_$unphased\_$variants/star_$parental\_withannotation/indexing.err
	
	echo "#!/bin/bash" > $codedir/MappingOptimization_$line\_runcommands.sh
	echo "" >> $codedir/MappingOptimization_$line\_runcommands.sh
	
	echo "# SLURM options" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --account=pfranken_bxd_map" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --job-name=mappingRun_"$line >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --time 0-03:00:00" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --nodes=1" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --ntasks=2" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --cpus-per-task=1" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --array=2-161%20" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --output=%x_%A-%a.out" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --error=%x_%A-%a.err" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --mem=27G" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --export=NONE" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --chdir "$mapdir >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --mail-user=nastassia.gobet@unil.ch" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "#SBATCH --mail-type=ALL" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "" >> $codedir/MappingOptimization_$line\_runcommands.sh
	
	echo "# load modules" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "module load HPC/Software" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "module load UHTS/Analysis/HTSeq/0.9.1" >> $codedir/MappingOptimization_$line\_runcommands.sh
	##echo "module load UHTS/Aligner/STAR/2.7.3a" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "" >> $codedir/MappingOptimization_$line\_runcommands.sh
	
	echo "# commands" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo 'date' >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "IN=\$(sed -n \${SLURM_ARRAY_TASK_ID}p $codedir/MappingOptimization_$line\_listcommands.sh)" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "echo Running analysis on row \${SLURM_ARRAY_TASK_ID}" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo "echo \$IN" >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo 'eval "$IN"' >> $codedir/MappingOptimization_$line\_runcommands.sh
	echo 'date' >> $codedir/MappingOptimization_$line\_runcommands.sh
	
	# copy reference directory from Franken NAS to Wally (!command to launch from 3009 computer and ask for login info!)
	##scp -r /mnt/nas/BXD/data/PersonalizedReferences/$line\_$unphased\_$variants ngobet@wally-front1.unil.ch:$refdir
	
	# run jobs
	##sbatch $codedir/MappingOptimization_$line\_runcommands.sh
	
	# remove reference directory copy on scratch
	##rm -r $refdir/$line\_$unphased\_$variants
done

# normalization [locally]

# data transfert on cluster
##############################################################################

# Part 2: eQTL analysis

# 20200622

# fastQTL is used to:
# tests association for a molecular phenotype only for variants that are 2 Mb above or below the transcription start site of the gene coding for the phenotype (--window 2e6, default is 1e6)
# chose to use seed one to help reproducibility of permutations (--seed 1, no default)
# uses the beta distribution approximation estimated from 1000 permutations to calculate adjusted p-values (--permute 1000, no default)
# EPACTS is used for compression and indexing of input files.

# load softwares
module load HPC/Software
module add UHTS/Analysis/FastQTL/2.184
module add UHTS/Analysis/EPACTS/3.2.6

# copy gene position file from 3009 (run command from Wally)
scp -r ngobet@pccig3009.unil.ch:/mnt/nas/BXD/references/transcriptome/GenePosition.txt $refdir/transcriptome_annotation/GenePosition.txt

# pre-processing: format and index genotypes
# transform genotypes to vcf format
scripts/transformGenotypesToVcf.py $refdir/genotypes/BXDGenotypes.geno $refdir/genotypes/BXDGenotypes.vcf
# compress and index genotypes
bgzip -f $refdir/genotypes/BXDGenotypes.vcf && tabix -p vcf $refdir/genotypes/BXDGenotypes.vcf.gz

# copy expression files on cluster from 3009 (command to run from cluster)
listsettings=($(ls -d $mapdir/* | cut -d "/" -f 11))
for setting in "${listsettings[@]}"
do
	echo $setting
	scp -r ngobet@pccig3009.unil.ch:/mnt/nas/BXD/data/transcriptome/2_mapping_STAR/MappingParametersOptimization/$setting/TMMnormalized_log2*PM_*.tab $mapdir/$setting/
done

listsettings=($(ls -d $mapdir/* | cut -d "/" -f 11))
for setting in "${listsettings[@]}"
do
	echo $setting
	mkdir $eQTLdir/$setting
	for tissue in Cortex Liver
	do
		echo $tissue
		for condition in NSD SD
		do
			echo $condition
			for $norm in CPM TPM
			do
				# transform phenotypes (gene expression) to bed format (UCSC format)
				scripts/transformGenePhenotypesToBedUCSC.py $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.tab $refdir/genotypes/BXDGenotypes.geno $refdir/transcriptome_annotation/GenePosition.txt $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.bed > $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition\_transformGenePhenotypesToBedUCSC.out 2> $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition\_transformGenePhenotypesToBedUCSC.err
				# compress and index phenotypes
				bgzip -f $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.bed && tabix -p bed $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.bed.gz
				
				# prepare commands
				fastQTL --vcf $refdir/genotypes/BXDGenotypes.vcf.gz --bed $mapdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.bed.gz --out $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition --commands 25 $codedir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.sh --window 2e6 --permute 1000 --seed 1
				# run commands
				bash $codedir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.sh > $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.out 2> $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.err
				# check if it run correctly, if not exclude last phenotype before error
				
				testvar=$(wc -l $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.err | cut -d " " -f 1)
				while [ $testvar -gt 0 ]
				do
					# check if it run correctly, if not exclude last phenotype before error
					echo "Error detected"
					# identify last phenotype
					issuecommandfull=$(grep "fastQTL" $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.err | sed 's/  \+/,/g' | cut -d "," -f 2)
					issuecommand=$(sed 's/ --exclude-phenotypes PhenoToExclude.txt//g' <<< $issuecommandfull)
					$issuecommand &> issueresults.txt
					grep "Processing gene" issueresults.txt | tail -n 1 | cut -d " " -f 3 | sed 's/[][]//g' >> $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition\_PhenoToExclude.txt
					# add exclude phenotype(s) option to script if not already contained
					if ! grep -q "\-\-exclude\-phenotypes" $codedir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.sh 
					then
						sed -i "s,$, --exclude-phenotypes $eQTLdir\/$setting\/TMMnormalized_log2$norm\_$tissue\_$condition\_PhenoToExclude.txt,g" $codedir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.sh
					fi
					# re-run excluding this phenotype
					bash $codedir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.sh > $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.out 2> $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.err
					rm issueresults.txt
					echo $testvar
					testvar=$(wc -l $eQTLdir/eQTL_$setting\_$tissue\_$condition\_$norm\_listcommands.err | cut -d " " -f 1)
				done
				echo "The eQTL detection successfully run."

				# post-processing
				# group results from different regions
				cat $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.chr* > $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.txt
				# add entry with NA for each gene excluded from analysis if any
				if [ -e $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition\_PhenoToExclude.txt ]
				then
					sed 's,$, 0 NA NA NA NA NA NA NA NA NA,g' $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition\_PhenoToExclude.txt >> $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.txt
				fi
				# clean up unneeded files
				rm $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition.chr* $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition
				# calculate q-values
				Rscript $codedir/correctMultiPhenotypeseQTL.R $eQTLdir/$setting/TMMnormalized_log2$norm\_$tissue\_$condition
			done
		done
	done
done

# copy expression files on 3009 (command to run from cluster)
scp -r $eQTLdir ngobet@pccig3009.unil.ch:/mnt/nas/BXD/data/transcriptome/

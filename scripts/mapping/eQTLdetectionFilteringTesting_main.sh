#!/bin/bash

###############################################################################

# CONTENT header

# eQTL detection for filtering evaluation

# GOAL: run cis eQTL detection with FastQTL for filtering evaluation

# DATE: 30 July 2020

# INPUT: normalized gene counts files (.tab)

# OUTPUT: eQTL statistics files (.txt)

###############################################################################

# directories definitions
refdir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/references
mapdir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/FilteringTesting/mapping
eQTLdir=/scratch/wally/FAC/FBM/CIG/pfranken/bxd_map2/data/FilteringTesting/eQTL
codedir=/users/ngobet/scripts

# transfert input on cluster
mkdir -p bxdmap_scratchdirectory2/data/FilteringTesting
mkdir -p $mapdir
mkdir -p $eQTLdir
scp -r ngobet@pccig3009.unil.ch:/mnt/nas/BXD/data/transcriptome/4_normalization/FilteringTesting/*.tab $mapdir

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
##scp -r ngobet@pccig3009.unil.ch:/mnt/nas/BXD/references/transcriptome/GenePosition.txt $refdir/transcriptome_annotation/GenePosition.txt

# pre-processing: format and index genotypes
# transform genotypes to vcf format
scripts/transformGenotypesToVcf.py $refdir/genotypes/BXDGenotypes.geno $refdir/genotypes/BXDGenotypes.vcf
# compress and index genotypes
bgzip -f $refdir/genotypes/BXDGenotypes.vcf && tabix -p vcf $refdir/genotypes/BXDGenotypes.vcf.gz

listexpressionfiles=($(ls -d $mapdir/*.tab | cut -d "/" -f 12))
for file in "${listexpressionfiles[@]}"
do
	echo $file
	filtering=$(basename -s .tab $file)
	echo $filtering
			# transform phenotypes (gene expression) to bed format (UCSC format)
			scripts/transformGenePhenotypesToBedUCSC.py $mapdir/$file $refdir/genotypes/BXDGenotypes.geno $refdir/transcriptome_annotation/GenePosition.txt $mapdir/$filtering.bed > $mapdir/$filtering\_transformGenePhenotypesToBedUCSC.out 2> $mapdir/$filtering\_transformGenePhenotypesToBedUCSC.err
			# compress and index phenotypes
			bgzip -f $mapdir/$filtering.bed && tabix -p bed $mapdir/$filtering.bed.gz
			
			# prepare commands
			fastQTL --vcf $refdir/genotypes/BXDGenotypes.vcf.gz --bed $mapdir/$filtering.bed.gz --out $eQTLdir/$filtering --commands 25 $codedir/eQTL_$filtering\_listcommands.sh --window 2e6 --permute 1000 --seed 1
			# run commands
			bash $codedir/eQTL_$filtering\_listcommands.sh > $eQTLdir/eQTL_$filtering\_listcommands.out 2> $eQTLdir/eQTL_$filtering\_listcommands.err
			
			# check if it run correctly, if not exclude last phenotype before error
			testvar=$(wc -l $eQTLdir/eQTL_$filtering\_listcommands.err | cut -d " " -f 1)
			while [ $testvar -gt 0 ]
			do
				# check if it run correctly, if not exclude last phenotype before error
				echo "Error detected"
				# identify last phenotype
				issuecommandfull=$(grep -m 1 "fastQTL" $eQTLdir/eQTL_$filtering\_listcommands.err | sed 's/  \+/,/g' | cut -d "," -f 2)
				issuecommand=$(sed 's/ --exclude-phenotypes PhenoToExclude.txt//g' <<< $issuecommandfull)
				$issuecommand &> issueresults.txt
				grep "Processing gene" issueresults.txt | tail -n 1 | cut -d " " -f 3 | sed 's/[][]//g' >> $eQTLdir/$filtering\_PhenoToExclude.txt
				# add exclude phenotype(s) option to script if not already contained
				if ! grep -q "\-\-exclude\-phenotypes" $codedir/eQTL_$filtering\_listcommands.sh 
				then
					sed -i "s,$, --exclude-phenotypes $eQTLdir\/$filtering\_PhenoToExclude.txt,g" $codedir/eQTL_$filtering\_listcommands.sh
				fi
				# re-run excluding this phenotype
				bash $codedir/eQTL_$filtering\_listcommands.sh > $eQTLdir/eQTL_$filtering\_listcommands.out 2> $eQTLdir/eQTL_$filtering\_listcommands.err
				rm issueresults.txt
				echo $testvar
				testvar=$(wc -l $eQTLdir/eQTL_$filtering\_listcommands.err | cut -d " " -f 1)
			done
			echo "The eQTL detection successfully run."

			# post-processing
			# group results from different regions
			cat $eQTLdir/$filtering.chr* > $eQTLdir/$filtering.txt
			# add entry with NA for each gene excluded from analysis if any
			if [ -e $eQTLdir/$filtering\_PhenoToExclude.txt ]
			then
				sed 's,$, 0 NA NA NA NA NA NA NA NA NA,g' $eQTLdir/$filtering\_PhenoToExclude.txt >> $eQTLdir/$filtering.txt
			fi
			# clean up unneeded files
			rm $eQTLdir/$filtering.chr* $eQTLdir/$filtering
			# calculate q-values
			Rscript $codedir/correctMultiPhenotypeseQTL.R $eQTLdir/$filtering
done

# copy files back on 3009 (command to run from cluster)
scp -r bxdmap_scratchdirectory2/data/FilteringTesting ngobet@pccig3009.unil.ch:/mnt/nas/BXD/data/transcriptome/

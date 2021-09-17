################
# Merge counts #
################

# GOAL: Merge gene expression counts from the different samples (one file per sample) to one file.

# DATE: 20200615

# AUTHOR: Nastassia Gobet

# INPUT
# 1) directory with samples counts. Example: data/transcriptome/2_mapping_STAR/OnPersonalizedRef/DBA_2J_nonrandomized_indelsSNVs/
# 2) suffix pattern to recognize samples counts files. Example: _permissive_ReadsPerGene.out.tab

# OUTPUT
# grouped counts file

#################################3

# CODE

# Variables definition
directory=$1
pattern=$2

# Control tests
# test if input samples counts files have the same number of lines (aka genes).
if [ -f $directory\Summary$pattern ]
then
	echo "Error: output file already exist: $directory\Summary$pattern"
	exit
fi

# Merge samples genes counts files
# get list of GeneID present in all files
totalsamples=$(ls $directory*$pattern | wc -l)
cut -f 1 $directory*$pattern > allids
sort allids | uniq -c | grep -w $totalsamples | sed 's/  */ /g' | cut -d " " -f 3 > consensusids
# creates a header file
echo -n "GeneID" > headers
# for each sample remove the gene ids (first column) to keep only the counts (second column) and save sample name in header
for samplecountfile in $directory*$pattern
do
	# Select number of column with counts
	col=$(awk '{print NF; exit}' $samplecountfile)
	echo "Taking counts from column $col (1-based)."
	# retrieve sample name
	sample=$(basename -s $pattern $samplecountfile)
	echo -ne "\t$sample" >> headers
	# retrieve geneid
	grep -f consensusids $samplecountfile | cut -f 1 > consensusidsordered
	# retrieve counts
	grep -f consensusids $samplecountfile | cut -f $col > values_$sample
done
# group all values in one file (first column is gene ids)
paste consensusidsordered values_* > body
# adds headers and remove the lines starting by __ or N_
cat headers <(echo) body | grep -v "^__" | grep -v "^N_" > $directory\Summary$pattern # <(echo) adds a end of line between the headers and the body
# remove temporary files
rm headers body values_* allids consensusids consensusidsordered

# Example of running this script: analysis/scripts/mapping/MergeCount.sh data/transcriptome/2_mapping_STAR/MappingParametersOptimization/genotypesandimputed_withoutannotation_Local_1_3/ _ReadsPerGene.out.tab

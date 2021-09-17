# GOAL: describe the proportion of DM genes in the different genes subtypes. (One of Olivier suggestion in my midthesis exam.)

# 7 January 2019

# 1) Define what are the different genes categories.
# 2) Count how many genes in each categories in total.
# 3) Represent results graphically.

# Needed files:
# annotation reference for B6 (mm10): references/transcriptome/Mus_musculus.GRCm38.94.gtf
# annotation reference for D2: references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.94.gtf
# list DM genes: data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_*.txt (Cortex, Liver)
# list of total genes used for DM analysis (genes with common gene names between B6 and D2): data/transcriptome/3_counting_htseq/OnGenome/commongenes.txt
#Beware this last file must have unix end of lines otherwise grep -f does not give any output.

# test

# list categories of genes
grep -w "gene" references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.94.gtf | grep -o "gene_biotype \"[a-zA-Z_]\{1,\}\"" | cut -d \" -f 2 | sort | uniq
grep -w "gene" references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep -o "gene_biotype \"[a-zA-Z_]\{1,\}\"" | cut -d \" -f 2 | sort | uniq

# count number of genes in each category (uniq -c)
grep -w "gene" references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.94.gtf | grep -o "gene_biotype \"[a-zA-Z_]\{1,\}\"" | cut -d \" -f 2 | sort | uniq -c
grep -w "gene" references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep -o "gene_biotype \"[a-zA-Z_]\{1,\}\"" | cut -d \" -f 2 | sort | uniq -c

# count for total genes with gene names
cut -d " " -f 2 convertD2 | sort | uniq -c 
cut -d " " -f 2 convertB6 | sort | uniq -c

###########

mkdir data/transcriptome/6_Biotypes

# Create a table with name to biotype category (to convert) for all genes with biotype and gene name.
grep -w "gene" references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.94.gtf | grep "gene_name" | cut -d \; -f 3,5 | cut -d \" -f 2,4 --output-delimiter " " > data/transcriptome/6_Biotypes/convertD2.txt
grep -w "gene" references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep "gene_name" | cut -d \; -f 3,5 | cut -d \" -f 2,4 --output-delimiter " " > data/transcriptome/6_Biotypes/convertB6.txt

# parse DM results to keep only DM genes (adjusted p-value < 0.05) and keep only gene names
awk '$6 < 0.05' data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_Cortex.txt | cut -f 1 > data/transcriptome/5_DifferentialMapping/DMgenenames_Cortex.txt
awk '$6 < 0.05' data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_Liver.txt | cut -f 1 > data/transcriptome/5_DifferentialMapping/DMgenenames_Liver.txt

# parse DM results to get all gene names used for DM analysis
cut -f 1 data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_Cortex.txt | sed '1d' > data/transcriptome/5_DifferentialMapping/genenames_Cortex.txt
cut -f 1 data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_Liver.txt | sed '1d' > data/transcriptome/5_DifferentialMapping/genenames_Liver.txt

# count genes in each biotype category for DM genes in Cortex
grep -w -f data/transcriptome/5_DifferentialMapping/DMgenenames_Cortex.txt data/transcriptome/6_Biotypes/convertB6.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryB6BiotypesDM_Cortex.tab
# count genes in each biotype category for DM genes in Liver
grep -w -f data/transcriptome/5_DifferentialMapping/DMgenenames_Liver.txt data/transcriptome/6_Biotypes/convertB6.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryB6BiotypesDM_Liver.tab
# count genes in each biotype category for genes considered for the DM analysis
grep -w -f data/transcriptome/5_DifferentialMapping/genenames_Cortex.txt data/transcriptome/6_Biotypes/convertB6.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryB6BiotypesTotal_Cortex.tab
grep -w -f data/transcriptome/5_DifferentialMapping/genenames_Liver.txt data/transcriptome/6_Biotypes/convertB6.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryB6BiotypesTotal_Liver.tab

# total genes all subcatergories
wc -l data/transcriptome/5_DifferentialMapping/DMgenenames_Cortex.txt data/transcriptome/5_DifferentialMapping/DMgenenames_Liver.txt data/transcriptome/5_DifferentialMapping/genenames_Cortex.txt data/transcriptome/5_DifferentialMapping/genenames_Liver.txt
# 4119 2497 13860 11132

# 24 January 2020

# GOAL: re-run biotypes analysis with D2 annotation

# count genes in each biotype category for DM genes in Cortex
grep -w -f data/transcriptome/5_DifferentialMapping/DMgenenames_Cortex.txt data/transcriptome/6_Biotypes/convertD2.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryD2BiotypesDM_Cortex.tab
# count genes in each biotype category for DM genes in Liver
grep -w -f data/transcriptome/5_DifferentialMapping/DMgenenames_Liver.txt data/transcriptome/6_Biotypes/convertD2.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryD2BiotypesDM_Liver.tab
# count genes in each biotype category for genes considered for the DM analysis
grep -w -f data/transcriptome/5_DifferentialMapping/genenames_Cortex.txt data/transcriptome/6_Biotypes/convertD2.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryD2BiotypesTotal_Cortex.tab
grep -w -f data/transcriptome/5_DifferentialMapping/genenames_Liver.txt data/transcriptome/6_Biotypes/convertD2.txt | cut -d " " -f 2 | sort | uniq -c > data/transcriptome/6_Biotypes/summaryD2BiotypesTotal_Liver.tab

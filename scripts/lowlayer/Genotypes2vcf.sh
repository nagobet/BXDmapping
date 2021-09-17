#!/bin/bash

#GOAL: Make a vcf file for all BXD + parental + F1 from mm10 genotypes on GeneNetwork

# 5 March 2020

# vcf files containing variant info
# SNPs: /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf
# indels: /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf
# new: ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/
# gvf file with some info on SVs (don't think i will be able to use them)
# SVs: /mnt/nas/BXD/data/genome/D2specificVariants/estd118_Keane_et_al_2011.2013-01-07.MGSCv37.gvf

# mm10 genotypes from GeneNetwork: /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGN.txt
# remove info lines at the beginning (starting by #,@ or ")
grep -v "^[#|@|\"]" /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGN.txt > geno_tmp
# remove columns with positions in CM or Mb
cut -f1-2,6- geno_tmp > geno_tmp2
# change H into 0/1, B into 0/0, D into 1/1, U into ./. (\b is the boundary of the word)
# help from: https://www.tutorialspoint.com/unix/unix-regular-expressions.htm and https://stackoverflow.com/questions/1032023/sed-whole-word-search-and-replace
sed -e 's/\bH\b/0\/1/g' -e 's/\bB\b/0\/0/g' -e 's/\bD\b/1\/1/g' -e 's/\bU\b/.\/./g' geno_tmp2 > geno_tmp3

# retrieve only ID of variants
cut -f 2 geno_tmp3 | grep -v "Locus" > IDvariants_tmp

# search variant ID from genotypes to vcf
grep -wf IDvariants_tmp /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > vcf_selected_tmp
# It takes a while (34 minutes).
grep -wf IDvariants_tmp /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf > vcfindels_selected_tmp

# number of variants ID
wc -l IDvariants_tmp
#7324
# number of SNP variants ID found
wc -l vcf_selected_tmp
#6852
# number of indels variants ID found
wc -l vcfindels_selected_tmp
#
# number of variants ID
cat vcf_selected_tmp vcfindels_selected_tmp | wc -l
#

### problem: when variant was not found, nothing is printed, so the lines don't match between the 2 files. => reselects matched variants in the genotype file (probably not the cleaner way of processing)
# retrieve only ID of variants that match between the genotypes and vcf
cut -f 3 vcf_selected_tmp > IDvariantsmatched_tmp
# reselects matched variants in the genotype file
grep -wf IDvariantsmatched_tmp geno_tmp3 > geno_tmp4

# remove the 2 first lines
cut -f3- geno_tmp4 > geno_tmp5

# remove the 10 column
cut -f1-9 vcf_selected_tmp > vcf_selected_tmp2
# replace the 9 th column
sed -e 's/GT:.\{1,\}/GT/g' vcf_selected_tmp2 > vcf_selected_tmp3

# retrieve header from vcf
grep "#CHROM" /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf | cut -f1-9 > header_vcf_tmp
# add header to vcf
cat header_vcf_tmp vcf_selected_tmp3> vcf_selected_tmp4

# retrieve header from genotypes
grep "Locus" geno_tmp | cut -f-5 > header_geno_tmp
# add header to genotypes
cat header_geno_tmp geno_tmp5 > geno_tmp6

# group files
paste vcf_selected_tmp4 geno_tmp6 > /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf

# remove temporary files
rm *tmp*

# replace / in strain names by _
sed -e 's/C57BL\/6J/C57BL_6J/g' -e 's/DBA\/2J/DBA_2J/g' -i /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf

##################

# 27 March 2020

# GOAL: trying to retrieve info on more variant ids from total SNPs and indels for all mouse strains.

# working directory
cd /mnt/nas/BXD

# mm10 genotypes from GeneNetwork: /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGN.txt
# remove info lines at the beginning (starting by #,@ or ")
grep -v "^[#|@|\"]" /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGN.txt > geno_tmp
# remove columns with positions in CM or Mb
cut -f1-2,6- geno_tmp > geno_tmp2
# change H into 0/1, B into 0/0, D into 1/1, U into ./. (\b is the boundary of the word)
# help from: https://www.tutorialspoint.com/unix/unix-regular-expressions.htm and https://stackoverflow.com/questions/1032023/sed-whole-word-search-and-replace
sed -e 's/\bH\b/0\/1/g' -e 's/\bB\b/0\/0/g' -e 's/\bD\b/1\/1/g' -e 's/\bU\b/.\/./g' geno_tmp2 > geno_tmp3

# retrieve only ID of variants
cut -f 2 geno_tmp3 | grep -v "Locus" > IDvariants_tmp

# search variant ID from genotypes to vcf
nohup bash -c 'grep -wf IDvariants_tmp data/genome/mgp.v6.merged.norm.snp.indels.sfiltered.vcf > vcf_selected_tmp' > nohupVCF.out 2> nohupVCF.err&

# number of variants ID
wc -l IDvariants_tmp
#7324
# number of variants ID found
wc -l vcf_selected_tmp
#6959

# find genotypes variantsID not found in vcf
sort IDvariants_tmp > tmpID
cut -f 3 vcf_selected_tmp | sort > tmpVCF
diff tmpID tmpVCF
diff tmpID tmpVCF | grep "<" | cut -d " " -f 2 > IDtofind
# number of variants in all strains file
grep -vc "#" data/genome/mgp.v6.merged.norm.snp.indels.sfiltered.vcf

######################################################################

# 8 April 2020

# try to retrieve more genotypes info with Maxime ideas
# 1) ftp://ftp.jax.org/dgatti/TPJ201600234
# 2) check if position corresponds to another variant.

# 1) search variant ID from genotypes to markers.txt
grep -wf IDtofind data/genome/markers.txt > markers_selected_tmp
wc -l markers_selected_tmp
#60

# 2) try to match by positions
grep "rs52712449" geno_tmp
grep "^15" data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > candidates

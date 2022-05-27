# Variants Characterization

#GOAL: characterize genomic variants D2-specific (vs genotypes from GeneNetwork)

# 30 March 2020

# Genotypes: /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGN.txt
# SNPs: /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf
# indels: /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf

# working directory
cd /mnt/nas/BXD

# how many SNVs are less than 100 bp from another SNV
grep -v "##" /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf | awk '{print>"/mnt/nas/BXD/data/genome/D2specificVariants/"$1".tab"}'
#help from: https://unix.stackexchange.com/questions/297683/splitting-a-file-into-multiple-files-based-on-1st-column-value

#run analysis/Mapping/SNVDistance.R

# remove temporary files
rm /mnt/nas/BXD/data/genome/D2specificVariants/*.tab

############################

# 1 April 2020

# annotation file: /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf
# Genotypes in vcf: /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf

# check for overlap between genotypes and genes
grep "^##" /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > info
cat info /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf > tmp.vcf
#Have to adding in vcf info of lines starting by ## to make it recognized by bedtools as a vcf file 
/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf -b tmp.vcf -wo > Genewithgenotypesoverlap&
/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf -b /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGNclean.vcf -u | grep -w gene > Genewithgenotypesoverlapmm10&
wc -l Genewithgenotypesoverlapmm10
#2738

# how many variants
grep -v -c "#" /mnt/nas/BXD/data/genome/BXD_Geno-19Jan2017_forGNclean.vcf
#7324


# how many variants genes have
###grep -w "gene" Genewithgenotypesoverlap | cut -f 10 | sort -g | uniq -c

# 22 April 2022

/home/ngobet/software/bedtools2/bin/intersectBed --help
# -u option if at least one everlap is found with B, reports once the entry in A
/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf -b /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf -u | grep -w gene > GenewithSNPoverlap&
wc -l GenewithSNPoverlap
#28930

/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf -b /mnt/nas/BXD/data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf -u | grep -w gene > Genewithindelsoverlap&
wc -l Genewithindelsoverlap
#23179

# how many are there in total genes
grep -w -c gene /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf
#54532


# how many genes overlap within D2 regions
/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf -b data/genome/BXD43_D2blocks.bed -u | grep -w gene > GenewithD2blocksoverlap&
wc -l GenewithD2blocksoverlap
#21457

### check overlap with exons???

# check overlap between mm9 genes and mm9 genotypes

grep "^#|@" data/genome/BXD_Geno_2001-2016.geno.txt > infomm9
cat info /mnt/nas/BXD/data/genome/GenotypesGNmm10.vcf > tmp.vcf

cat data/genome/BXD_Geno_2001-2016.geno.txt > tmpmm9.vcf
# check exons overlapping
/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/references/transcriptome/Mus_musculus.NCBIM37.67.gtf -b /mnt/nas/BXD/data/genome/BXD_Geno_2001-2016clean.vcf -u | grep -w exon > Exonwithgenotypesoverlapmm9&
# determine unique genes overlaped
cut -f 9 Exonwithgenotypesoverlapmm9 | cut -d ";" -f 1 | sort | uniq -c | wc -l
#555

# how many genotypes overlap
/home/ngobet/software/bedtools2/bin/intersectBed -a /mnt/nas/BXD/data/genome/BXD_Geno_2001-2016clean.vcf -b /mnt/nas/BXD/references/transcriptome/Mus_musculus.NCBIM37.67.gtf -u | sort | uniq | wc -l
#554

# total number of variants
grep -v -c "#" /mnt/nas/BXD/data/genome/BXD_Geno_2001-2016clean.vcf
#3811

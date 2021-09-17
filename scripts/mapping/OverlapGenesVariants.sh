# Overlap genes and variants

# Goal: find if genes with variants are more differentially mapped. If SNPs or indels have more effects.

# Date: 9 June 2021

# get list of expressed genes
cut -f 1 data/transcriptome/5_DifferentialMapping/DifferentialMappingBXDvsB6mm10permissive_Limma_Liver.txt | grep -v "logFC" > liver_expressedgenes
cut -f 1 data/transcriptome/5_DifferentialMapping/DifferentialMappingBXDvsB6mm10permissive_Limma_Cortex.txt | grep -v "logFC" > cortex_expressedgenes

# retrieve genes coordinates for "gene" entries
sed 's/.*/"&"/' liver_expressedgenes > liver_expressedgenes.txt
sed 's/.*/"&"/' cortex_expressedgenes > cortex_expressedgenes.txt
nohup bash -c 'grep -w gene references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep -w -f liver_expressedgenes.txt > liver.gtf' > liver_nohup.out 2> liver_nohup.err&
nohup bash -c 'grep -w gene references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep -w -f cortex_expressedgenes.txt > cortex.gtf' > cortex_nohup.out 2> cortex_nohup.err&
### or better exons only?

# remove heterozygous variants from list of variants
grep -v -w "Het" data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > snps.vcf
grep -v -w "Het" data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf > indels.vcf

# retrieve overlap variants per gene
~/software/bedtools2/bin/intersectBed -c -a liver.gtf -b snps.vcf > data/transcriptome/5_DifferentialMapping/liver_genesSNVsoverlap_BXDvsB6mm10permissive.txt&
~/software/bedtools2/bin/intersectBed -c -a liver.gtf -b indels.vcf > data/transcriptome/5_DifferentialMapping/liver_genesindelsoverlap_BXDvsB6mm10permissive.txt&
~/software/bedtools2/bin/intersectBed -c -a cortex.gtf -b snps.vcf > data/transcriptome/5_DifferentialMapping/cortex_genesSNVsoverlap_BXDvsB6mm10permissive.txt&
~/software/bedtools2/bin/intersectBed -c -a cortex.gtf -b indels.vcf > data/transcriptome/5_DifferentialMapping/cortex_genesindelsoverlap_BXDvsB6mm10permissive.txt&

# clean temporary files
rm liver_expressedgenes* cortex_expressedgenes* liver.gtf cortex.gtf *_nohup.* snps.vcf indels.vcf

###############################################

# do the same for comparing genes with D2 assembly and GRCm38

# get list of expressed genes
cut -f 1 data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_Liver_20190917.txt | grep -v "logFC" > liver_expressedgenes
cut -f 1 data/transcriptome/5_DifferentialMapping/DifferentialMapping_Limma_Cortex_20190917.txt | grep -v "logFC" > cortex_expressedgenes

# retrieve genes coordinates for "gene" entries
sed 's/.*/"&"/' liver_expressedgenes > liver_expressedgenes.txt
sed 's/.*/"&"/' cortex_expressedgenes > cortex_expressedgenes.txt
nohup bash -c 'grep -w gene references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep -w -f liver_expressedgenes.txt > liver.gtf' > liver_nohup.out 2> liver_nohup.err&
nohup bash -c 'grep -w gene references/transcriptome/Mus_musculus.GRCm38.94.gtf | grep -w -f cortex_expressedgenes.txt > cortex.gtf' > cortex_nohup.out 2> cortex_nohup.err&
### or better exons only?

# remove heterozygous variants from list of variants
grep -v -w "Het" data/genome/D2specificVariants/DBA_2J.mgp.v5.snps.dbSNP142.vcf > snps.vcf
grep -v -w "Het" data/genome/D2specificVariants/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf > indels.vcf

# retrieve overlap variants per gene
~/software/bedtools2/bin/intersectBed -c -a liver.gtf -b snps.vcf > data/transcriptome/5_DifferentialMapping/liver_genesSNVsoverlap_D2assemblyvsB6mm10.txt&
~/software/bedtools2/bin/intersectBed -c -a liver.gtf -b indels.vcf > data/transcriptome/5_DifferentialMapping/liver_genesindelsoverlap_D2assemblyvsB6mm10.txt&
~/software/bedtools2/bin/intersectBed -c -a cortex.gtf -b snps.vcf > data/transcriptome/5_DifferentialMapping/cortex_genesSNVsoverlap_D2assemblyvsB6mm10.txt&
~/software/bedtools2/bin/intersectBed -c -a cortex.gtf -b indels.vcf > data/transcriptome/5_DifferentialMapping/cortex_genesindelsoverlap_D2assemblyvsB6mm10.txt&

# clean temporary files
rm liver_expressedgenes* cortex_expressedgenes* liver.gtf cortex.gtf *_nohup.* snps.vcf indels.vcf

# GOAL: Mapping with Kallisto in local 

# generate index B6 mm10
##/software/Alignment/Kallisto/0.45.0/bin/kallisto index -i /mnt/nas/BXD/references/transcriptome/kallisto_mm10 /mnt/nas/BXD/references/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa
# generate index D2
##/software/Alignment/Kallisto/0.45.0/bin/kallisto index -i /mnt/nas/BXD/references/transcriptome/kallisto_D2 /mnt/nas/BXD/references/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa
#core dumps

# take the one that Maxime made
# path to index 
/index/Kallisto/mm10_ensembl_95.idx
# path to gtf
/reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf

# command quantification
/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /index/Kallisto/mm10_ensembl_95.idx -o /home/ngobet/projects/BXD/data/TMP --single -l 200 -s 200 --gtf /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/05nsd.fastq.gz &
/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /index/Kallisto/mm10_ensembl_95.idx -o /home/ngobet/projects/BXD/data/TMP/L05nsd --single -l 200 -s 200 --gtf /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L05nsd.fastq.gz &

/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /index/Kallisto/mm10_ensembl_95.idx -o /home/ngobet/projects/BXD/data/TMP/L89sd --single -l 200 -s 200 --gtf /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/L89sd.fastq.gz &
/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /index/Kallisto/mm10_ensembl_95.idx -o /home/ngobet/projects/BXD/data/TMP/89sd --single -l 200 -s 200 --gtf /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/89sd.fastq.gz &
/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /index/Kallisto/mm10_ensembl_95.idx -o /home/ngobet/projects/BXD/data/TMP/89sd --single -l 200 -s 20 --gtf /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/filtered/89sd.fastq.gz &

# re- build from Maxime fasta files
/software/Alignment/Kallisto/0.45.0/bin/kallisto index -i mm10_ensembl_95.idx /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.cdna.all.fa.gz &
#working

# launch the jobs to quantify with Kallisto
for sample in $(<data/transcriptome/RNAseq/listSamplesGSM.txt)
do
  ##name=$(basename $list)
  #echo $name
  #cat $list
  /software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /index/Kallisto/mm10_ensembl_95.idx -o data/transcriptome/2.5_quantification_Kallisto/$sample/ --single -l 200 -s 30 --gtf /reference/transcriptome/mmusculus/mm10/ensembl/release-95/Mus_musculus.GRCm38.95.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/GEO/$sample.fastq.gz
done

# retrieve mapping stats
grep "n_pseudoaligned" data/transcriptome/2.5_quantification_Kallisto/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/n_pseudoaligned.tab
grep "n_unique" data/transcriptome/2.5_quantification_Kallisto/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/n_unique.tab
grep "n_processed" data/transcriptome/2.5_quantification_Kallisto/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/n_processed.tab

# 12 December 2019

# build my own index (to have the same reference that I used for my tests with STAR.
mkdir /home/ngobet/projects/BXD/references
# generate index D2
##nohup /software/Alignment/Kallisto/0.45.0/bin/kallisto index -i /home/ngobet/projects/BXD/references/kallisto_D2.idx /mnt/nas/BXD/references/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa.gz &
# generate index B6 mm10
##/software/Alignment/Kallisto/0.45.0/bin/kallisto index -i /home/ngobet/projects/BXD/references/kallisto_B6mm10.idx /mnt/nas/BXD/references/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz &
# try to compress the files before to run indexing
##gzip -k /mnt/nas/BXD/references/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa
##gzip -k /mnt/nas/BXD/references/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa &
# core dumps (does not change if .fa or .fa.gz) => I think I have to give the fasta of transcripts only (not full DNA fasta)
# "To build the indices, download the full transcriptomes from Ensembl (files ending in cdna.all.fa.gz) and build the indices with kallisto." (source: https://github.com/pachterlab/kallisto-transcriptome-indices)
# download from Ensembl
# ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
# ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus_dba2j/cdna/Mus_musculus_dba2j.DBA_2J_v1.cdna.all.fa.gz
# generate index B6 mm10
nohup /software/Alignment/Kallisto/0.45.0/bin/kallisto index -i /mnt/nas/BXD/references/transcriptome/kallisto_B6mm10.idx /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.cdna.all.fa.gz > /mnt/nas/BXD/references/transcriptome/kallisto_B6mm10_idx_building_info.txt &
# generate index D2
nohup /software/Alignment/Kallisto/0.45.0/bin/kallisto index -i /mnt/nas/BXD/references/transcriptome/kallisto_D2.idx /mnt/nas/BXD/references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.cdna.all.fa.gz > /mnt/nas/BXD/references/transcriptome/kallisto_D2_idx_building_info.txt &

# launch the jobs to quantify with Kallisto
mkdir data/transcriptome/2.5_quantification_Kallisto/D2/
mkdir data/transcriptome/2.5_quantification_Kallisto/B6mm10/
for sample in $(<data/transcriptome/RNAseq/listSamplesGSM.txt)
do
	/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /mnt/nas/BXD/references/transcriptome/kallisto_D2.idx -o data/transcriptome/2.5_quantification_Kallisto/D2/$sample/ --single -l 200 -s 30 --gtf /mnt/nas/BXD/references/transcriptome/Mus_musculus_dba2j.DBA_2J_v1.94.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/GEO/$sample.fastq.gz
	/software/Alignment/Kallisto/0.45.0/bin/kallisto quant -i /mnt/nas/BXD/references/transcriptome/kallisto_B6mm10.idx -o data/transcriptome/2.5_quantification_Kallisto/B6mm10/$sample/ --single -l 200 -s 30 --gtf /mnt/nas/BXD/references/transcriptome/Mus_musculus.GRCm38.94.gtf --genomebam -t 10 /mnt/nas/BXD/data/transcriptome/RNAseq/GEO/$sample.fastq.gz
done

# 13 December 2019

# retrieve mapping stats
grep "n_pseudoaligned" data/transcriptome/2.5_quantification_Kallisto/D2/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/D2/n_pseudoaligned.tab
grep "n_unique" data/transcriptome/2.5_quantification_Kallisto/D2/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/D2/n_unique.tab
grep "n_processed" data/transcriptome/2.5_quantification_Kallisto/D2/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/D2/n_processed.tab

grep "n_pseudoaligned" data/transcriptome/2.5_quantification_Kallisto/B6mm10/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_pseudoaligned.tab
grep "n_unique" data/transcriptome/2.5_quantification_Kallisto/B6mm10/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_unique.tab
grep "n_processed" data/transcriptome/2.5_quantification_Kallisto/B6mm10/GSM*/*.json | cut -d : -f 3 | sed 's/[ ,]//g' > data/transcriptome/2.5_quantification_Kallisto/B6mm10/n_processed.tab

# get assembly statistics

# forum: https://www.biostars.org/p/393642/
# software: https://github.com/NBISweden/GAAS

# do inside an environnment as said here (https://exerror.com/found-conflicts-looking-for-incompatible-packages-this-can-take-several-minutes-press-ctrl-c-to-abort/)
conda create --name myenv
conda activate myenv
gaas_fasta_statistics.pl -f references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa > mm10_assembly_stats.txt&
gaas_fasta_statistics.pl -f references/genome/Mus_musculus_dba2j.DBA_2J_v1.dna_sm.toplevel.fa > D2_assembly_stats.txt&
# check for mm9 too
gaas_fasta_statistics.pl -f  references/genome/Mus_musculus.NCBIM37.67.dna.toplevel.fa > mm9_assembly_stats.txt&

gaas_fasta_statistics.pl -f references/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa > mm10_assembly_stats.txt&

# for mm10 I don't get the same results than Paul nor that I got last time. ??? Maybe not run on the same file? (not same patches, scaffolds, haplotypes included?)
gaas_fasta_statistics.pl -f references/genome/Mus_musculus.GRCm38.dna.toplevel.fa > mm10other_assembly_stats.txt&

# format in a table results
echo "" > tmp_emp
echo "GRCm38" > tmp_G
echo "D2" > tmp_D
paste tmp_emp tmp_G tmp_D > tmp_header  

cut -d "|" -f 2  mm10_assembly_stats.txt | grep -v "\\-\\-\\-" | grep -v "Ana" | grep -v "Mus" > tmp_stats
cut -d "|" -f 3  mm10_assembly_stats.txt | grep -v "\\-\\-\\-" | grep . > tmp_mm10
cut -d "|" -f 3  D2_assembly_stats.txt | grep -v "\\-\\-\\-" | grep . > tmp_D2
##wc -l tmp_*

paste tmp_stats tmp_mm10 tmp_D2 > tmp_body

cat tmp_header tmp_body > references/genome/AssemblyStatistics.tsv

rm tmp_*

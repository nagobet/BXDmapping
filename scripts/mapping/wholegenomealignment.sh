# whole genome alignment

# mummer
#https://mummer4.github.io/tutorial/tutorial.html


mkdir data/genome/wholegenomealignment
/home/ngobet/software/mummer-3.9.4alpha/bin/mummer
/home/ngobet/software/mummer-3.9.4alpha/bin/mummer -mum -b -c references/Mus_musculus.GRCm38.dna.major_chromosomes.fa references/Mus_musculus_dba2j.DBA_2J_v1.dna.major_chromosomes.fa > data/genome/wholegenomealignment/mummer.mums
### not working

# Blast by chromosome
# 20190806

# select chromosome 1
grep -n ">" Mus_musculus_dba2j.DBA_2J_v1.dna.major_chromosomes.fa
#1:>1 dna:chromosome chromosome:DBA_2J_v1:1:1:192928109:1 REF
#3215471:>2 dna:chromosome chromosome:DBA_2J_v1:2:1:179964709:1 REF
#6214884:>3 dna:chromosome chromosome:DBA_2J_v1:3:1:157190453:1 REF
8834726:>4 dna:chromosome chromosome:DBA_2J_v1:4:1:153890123:1 REF
11399563:>5 dna:chromosome chromosome:DBA_2J_v1:5:1:150592947:1 REF
13909447:>6 dna:chromosome chromosome:DBA_2J_v1:6:1:147407711:1 REF
16366244:>7 dna:chromosome chromosome:DBA_2J_v1:7:1:144344007:1 REF
18771979:>8 dna:chromosome chromosome:DBA_2J_v1:8:1:126593561:1 REF
20881873:>9 dna:chromosome chromosome:DBA_2J_v1:9:1:122266228:1 REF
22919645:>10 dna:chromosome chromosome:DBA_2J_v1:10:1:128213663:1 REF
25056541:>11 dna:chromosome chromosome:DBA_2J_v1:11:1:120239594:1 REF
27060536:>12 dna:chromosome chromosome:DBA_2J_v1:12:1:117921364:1 REF
29025894:>13 dna:chromosome chromosome:DBA_2J_v1:13:1:117427490:1 REF
30983020:>14 dna:chromosome chromosome:DBA_2J_v1:14:1:117845551:1 REF
32947114:>15 dna:chromosome chromosome:DBA_2J_v1:15:1:101173774:1 REF
34633345:>16 dna:chromosome chromosome:DBA_2J_v1:16:1:94787448:1 REF
36213137:>17 dna:chromosome chromosome:DBA_2J_v1:17:1:93736615:1 REF
37775415:>18 dna:chromosome chromosome:DBA_2J_v1:18:1:87488644:1 REF
39233561:>19 dna:chromosome chromosome:DBA_2J_v1:19:1:58689003:1 REF
40211713:>X dna:chromosome chromosome:DBA_2J_v1:X:1:166092115:1 REF


sed -n '1,3215470p;3215471q' /f/BXD/references/Mus_musculus_dba2j.DBA_2J_v1.dna.major_chromosomes.fa > /c/Users/n_gob/Desktop/D2_chr1.fa
sed -n '30983020,32947113p;32947114q' /f/BXD/references/Mus_musculus_dba2j.DBA_2J_v1.dna.major_chromosomes.fa > /c/Users/n_gob/Desktop/D2_chr14.fa

 
grep -n ">" Mus_musculus.GRCm38.dna.major_chromosomes.fa
1:>1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
3257869:>10 dna:chromosome chromosome:GRCm38:10:1:130694993:1 REF
5436120:>11 dna:chromosome chromosome:GRCm38:11:1:122082543:1 REF
7470831:>12 dna:chromosome chromosome:GRCm38:12:1:120129022:1 REF
9472983:>13 dna:chromosome chromosome:GRCm38:13:1:120421639:1 REF
11480012:>14 dna:chromosome chromosome:GRCm38:14:1:124902244:1 REF
13561718:>15 dna:chromosome chromosome:GRCm38:15:1:104043685:1 REF
15295781:>16 dna:chromosome chromosome:GRCm38:16:1:98207768:1 REF
16932579:>17 dna:chromosome chromosome:GRCm38:17:1:94987271:1 REF
18515702:>18 dna:chromosome chromosome:GRCm38:18:1:90702639:1 REF
20027414:>19 dna:chromosome chromosome:GRCm38:19:1:61431566:1 REF
21051275:>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF
24086497:>3 dna:chromosome chromosome:GRCm38:3:1:160039680:1 REF
26753826:>4 dna:chromosome chromosome:GRCm38:4:1:156508116:1 REF
29362296:>5 dna:chromosome chromosome:GRCm38:5:1:151834684:1 REF
31892876:>6 dna:chromosome chromosome:GRCm38:6:1:149736546:1 REF
34388487:>7 dna:chromosome chromosome:GRCm38:7:1:145441459:1 REF
36812513:>8 dna:chromosome chromosome:GRCm38:8:1:129401213:1 REF
38969201:>9 dna:chromosome chromosome:GRCm38:9:1:124595110:1 REF


sed -n '1,3257868p;3257869q' /f/BXD/references/Mus_musculus.GRCm38.dna.major_chromosomes.fa > /c/Users/n_gob/Desktop/B6_chr1.fa
sed -n '11480012,13561717p;13561718q' /f/BXD/references/Mus_musculus.GRCm38.dna.major_chromosomes.fa > /c/Users/n_gob/Desktop/B6_chr14.fa
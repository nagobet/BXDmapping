# 20 December 2019

# GOAL: erase bam files for mismatches 6, 7, 8, and 9 to avoid memory issues

# testing
##rm data/transcriptome/2_mapping_STAR/OnTranscriptome/D2/*[6789]mismatch*.bam data/transcriptome/2_mapping_STAR/OnTranscriptome/B6mm10/*[6789]mismatch*.bam
##watch -n 10 ls data/transcriptome/2_mapping_STAR/OnTranscriptome/D2/*[67891]mismatch*.bam data/transcriptome/2_mapping_STAR/OnTranscriptome/B6mm10/*[6789]mismatch*.bam >> tmp
##nohup watch -n 10 date >> tmp &
##nohup watch -n 30 rm *[6789]mismatch*.txt &

# final command
nohup watch -n 3600 rm data/transcriptome/2_mapping_STAR/OnTranscriptome/D2/*[6789]mismatch*.bam data/transcriptome/2_mapping_STAR/OnTranscriptome/B6mm10/*[6789]mismatch*.bam &
# to stop: kill 31361
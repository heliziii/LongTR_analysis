#PBS -q hotel 
#PBS -N merge
#PBS -l nodes=1:ppn=20:cascade
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group

cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont

 

#samtools merge -@ 20 HG002_ont_DUPLEX_hg38_all.bam HG002_ont_DUPLEX_hg38_sorted_{1..12}.bam

samtools sort -@ 20  HG002_ont_DUPLEX_hg38_all.bam -o HG002_ont_DUPLEX_hg38_sorted_all_2.bam

samtools index HG002_ont_DUPLEX_hg38_sorted_all_2.bam

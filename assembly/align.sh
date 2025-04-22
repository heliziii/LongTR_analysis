#PBS -q hotel 
#PBS -N HG002_align
#PBS -l nodes=1:ppn=5
#PBS -l walltime=04:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/sequence_compare

time /projects/ps-gymreklab/helia/packages/minimap2/minimap2 -ax asm5 /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa hg002v1.0.fasta -t 5 > HG002_assembly_hg38.sam  

samtools view -bh HG002_assembly_hg38.sam > HG002_assembly_hg38.bam

samtools sort HG002_assembly_hg38.bam -o HG002_assembly_hg38_sorted.bam

samtools index HG002_assembly_hg38_sorted.bam 

rm HG002_assembly_hg38.sam HG002_assembly_hg38.bam

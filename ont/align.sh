#PBS -q hotel 
#PBS -N align
#PBS -l nodes=1:ppn=20:cascade
#PBS -l walltime=100:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group

cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"

echo $index

wget https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/11_15_22_R1041_Duplex_HG002_"$index"_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz
 
time /projects/ps-gymreklab/helia/packages/minimap2/minimap2 -ax map-ont -t 20 $ref 11_15_22_R1041_Duplex_HG002_"$index"_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz > HG002_ont_DUPLEX_hg38_"$index".sam

samtools view -bh HG002_ont_DUPLEX_hg38_"$index".sam > HG002_ont_DUPLEX_hg38_"$index".bam

samtools sort HG002_ont_DUPLEX_hg38_"$index".bam -o HG002_ont_DUPLEX_hg38_sorted_"$index".bam

samtools index HG002_ont_DUPLEX_hg38_sorted_"$index".bam

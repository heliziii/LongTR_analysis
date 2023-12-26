#PBS -q hotel 
#PBS -N trgt_hg002
#PBS -l nodes=1:ppn=10:cascade
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group

cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"

reads=HG002_ont_DUPLEX_hg38_sorted_all.bam
sample=HG002

time /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/trgtv0.5.0/trgt-v0.5.0-linux_x86_64 --genome $ref \
                                --repeats /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/adotto_repeats.hg38.bed \
                                --reads $reads \
                                --output-prefix "$sample"_wg_trgt_hg38 \
                                --karyotype XY \
                                --flank-len 25 \
                                --max-depth 10000 \
				--min-read-quality 0.000000001 \
				--threads 10 

bcftools sort -Ob -o "$sample"_wg_trgt_hg38.sorted.vcf.gz "$sample"_wg_trgt_hg38.vcf.gz
bcftools index "$sample"_wg_trgt_hg38.sorted.vcf.gz

#PBS -q hotel 
#PBS -N trgt
#PBS -l nodes=1:ppn=12
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/homopolymer_fix/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/homopolymer_fix/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/homopolymer_fix
#reads="/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG003.m84010_220919_235306_s2.GRCh38.bam"
reads="/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG004.m84010_220919_232145_s1.GRCh38.bam"
#reads="/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG002.m84005_220827_014912_s1.GRCh38.bam"
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"
sample=HG004
time /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/trgtv0.5.0/trgt-v0.5.0-linux_x86_64 --genome $ref \
                                --repeats TRGT_homopolymer_ref.bed \
                                --reads $reads \
                                --output-prefix TRGT_calls/"$sample"_wg_trgt_hg38 \
                                --karyotype XY \
                                --flank-len 25 \
                                --max-depth 10000 \
				--threads 10

bcftools sort -Ob -o TRGT_calls/"$sample"_homo_trgt_hg38.sorted.vcf.gz TRGT_calls/"$sample"_wg_trgt_hg38.vcf.gz
bcftools index TRGT_calls/"$sample"_homo_trgt_hg38.sorted.vcf.gz


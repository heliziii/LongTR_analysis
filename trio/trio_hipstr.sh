#PBS -q hotel 
#PBS -N trio_hipstr
#PBS -l nodes=1:ppn=1:cascade
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/trio
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"

echo $chr
/projects/ps-gymreklab/helia/HipSTR_LR/LongSTR/LongTR \
--bams reads/HG002.m84005_220827_014912_s1.GRCh38.bam,reads/HG003.m84010_220919_235306_s2.GRCh38.bam,reads/HG004.m84010_220919_232145_s1.GRCh38.bam \
--fasta $ref \
--regions /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/hipstr_adotto_set/hipstr_catalog_adotto_chr"$chr".bed \
--tr-vcf calls/trio_hipstr_hg38_chr"$chr"_dec15.vcf.gz \
--bam-samps HG002,HG003,HG004 --bam-libs HG002,HG003,HG004 \
--min-reads 5 \
--output-filters \
--max-tr-len 10000 \
--log log/log_"$chr"_dec15.txt \
--min-sum-qual -1e18 \
--skip-assembly \
--phased-bam

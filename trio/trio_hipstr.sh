#!/bin/bash -l

#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J trio
#SBATCH -t 20:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/log/gt-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/log/gt-%j.err-%N
#SBATCH -N 1
#SBATCH -V
#SBATCH -n 1
#SBATCH --mem 4G
#SBATCH --export ALL

cd /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio
ref="/tscc/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"
chr=$1
echo $chr
/tscc/projects/ps-gymreklab/helia/HipSTR_LR/LongSTR/LongTR \
--bams reads/HG002.m84005_220827_014912_s1.GRCh38.bam,reads/HG003.m84010_220919_235306_s2.GRCh38.bam,reads/HG004.m84010_220919_232145_s1.GRCh38.bam \
--fasta $ref \
--regions /tscc/projects/ps-gymreklab/helia/HipSTR_LR/TRGT/hipstr_adotto_set/hipstr_catalog_adotto_chr"$chr".bed \
--tr-vcf calls/trio_hipstr_hg38_chr"$chr"_apr23.vcf.gz \
--bam-samps HG002,HG003,HG004 --bam-libs HG002,HG003,HG004 \
--min-reads 5 \
--output-filters \
--max-tr-len 10000 \
--log log/log_"$chr"_apr23.txt \
--skip-assembly \
--phased-bam

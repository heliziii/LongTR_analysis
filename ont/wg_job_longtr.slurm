#!/bin/bash


#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J gt
#SBATCH -t 30:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log/gt-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log/gt-%j.err-%N
#SBATCH -N 1
#SBATCH -V
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --export ALL

cd /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/ont
ref="/tscc/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"
chr="$1"

module load samtools

time /tscc/projects/ps-gymreklab/helia/HipSTR_LR/LongSTR/LongTR \
--bams simplex/hg002.haplotagged.chr21.bam \
--fasta $ref \
--regions /tscc/projects/ps-gymreklab/helia/HipSTR_LR/TRGT/hipstr_adotto_set/hipstr_catalog_adotto.bed \
--tr-vcf review_output/HG002_wg_longtr_ont_march_chr"$chr"_simplex.vcf.gz \
--bam-samps HG002 --bam-libs HG002 \
--min-reads 4 \
--chrom chr"$chr" \
--output-filters \
--haploid-chrs chrX,chrY \
--max-tr-len 10000 \
--min-mean-qual 20 \
--skip-assembly \
--phased-bam \
--indel-flank-len 25 \
--log log/longtr_march17_chr"$chr"_simplex.txt

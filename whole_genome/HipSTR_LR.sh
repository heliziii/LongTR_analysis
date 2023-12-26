#!/bin/bash

#PBS -q hotel 
#PBS -N HipSTR_LR
#PBS -l nodes=1:ppn=1:cascade
#PBS -l walltime=50:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group
#PBS -V


dir=/projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome
echo $chr

time /projects/ps-gymreklab/helia/HipSTR_LR/LongSTR/HipSTR \
--bams /projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG002.m84005_220827_014912_s1.GRCh38.bam \
--fasta /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
--regions "$dir"/metadata/GRCh38.hipstr_reference_pre.bed \
--tr-vcf "$dir"/output/HG002_wg_HiFi_HipSTR_long_chr"$chr".vcf.gz \
--bam-samps HG002 --bam-libs HG002 \
--min-reads 5 \
--output-filters \
--max-tr-len 10000 \
--min-sum-qual -1e18 \
--skip-assembly \
--haploid-chrs chrX,chrY \
--chrom chr"$chr" \
--phased-bam

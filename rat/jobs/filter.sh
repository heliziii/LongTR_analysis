#!/bin/bash

#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J dumpstr_rat
#SBATCH -t 2:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/log/dumpstr-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/log/dumpstr-%j.err-%N
#SBATCH -N 1
#SBATCH -V
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --export ALL

source /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/new_env/bin/activate
n=$1
address=/tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat
time dumpSTR  \
    --vcf $address/output/rat_"$n".vcf.gz \
    --out $address/filtered_calls/rat_"$n"_filtered \
    --vcftype hipstr \
    --min-locus-callrate 0.75 \
    --filter-regions  /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/parascopy_files/parascopy_segdups_filtered.bed.gz \
    --filter-regions-names SEGDUP \
    --drop-filtered &&

bgzip -f $address/filtered_calls/rat_"$n"_filtered.vcf && 

tabix -p vcf $address/filtered_calls/rat_"$n"_filtered.vcf.gz

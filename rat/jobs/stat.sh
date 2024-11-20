#!/bin/bash

#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J statstr_rat
#SBATCH -t 2:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/log/statstr-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/log/statstr-%j.err-%N
#SBATCH -N 1
#SBATCH -V
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --export ALL

cd /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/
source /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/new_env/bin/activate
n=$1
address=/tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat

bcftools filter -i 'INFO/PERIOD!=1' $address/filtered_calls/rat_"$n"_filtered.vcf.gz -o $address/filtered_calls/rat_"$n"_filtered_non_homo.vcf.gz -O z
tabix -p vcf $address/filtered_calls/rat_"$n"_filtered_non_homo.vcf.gz
time statSTR  \
    --vcf $address/filtered_calls/rat_"$n"_filtered_non_homo.vcf.gz \
    --afreq \
    --nalleles \
    --het \
    --mean \
    --vcftype hipstr \
    --out stats/stats_rat_"$n"_non_homo 

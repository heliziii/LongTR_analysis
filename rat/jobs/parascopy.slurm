#!/bin/bash
#SBATCH -q hotel
#SBATCH -p hotel
#SBATCH -J parascopy
#SBATCH -t 20:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/log/slurm-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/log/slurm-%j.err-%N
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 13

module load cpu/0.17.3
module load anaconda3/2021.05-xvr4x2l

source /cm/shared/apps/spack/0.17.3/cpu/opt/spack/linux-rocky9-cascadelake/gcc-11.2.0/anaconda3-2021.05-xvr4x2lssxb4duve5743jja5be44aale/bin/activate  paras_env

cd /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/

#parascopy pretable -f GCF_015227675.2_mRatBN7.2_genomic.fna -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/parascopy_files/pretable.bed.gz -@ 12 --tmp-dir /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/parascopy_files/temp
parascopy table -i /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/parascopy_files/pretable_sorted.bed.gz -f GCF_015227675.2_mRatBN7.2_genomic.fna -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/rat/parascopy_files/table.bed.gz -@ 12 


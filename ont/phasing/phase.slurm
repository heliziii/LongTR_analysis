#!/bin/bash


#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J phasing
#SBATCH -t 30:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log/phase-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log/phase-%j.err-%N
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5

dir=/tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/ont

source /cm/shared/apps/spack/0.17.3/cpu/opt/spack/linux-rocky9-cascadelake/gcc-11.2.0/anaconda3-2021.05-xvr4x2lssxb4duve5743jja5be44aale/bin/activate
conda activate whatshap-env

whatshap haplotag  --output "$dir"/phasing/HG002_ont_DUPLEX_hg38_sorted_haplotagged_all.bam \
	          --reference /tscc/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
	          "$dir"/phasing/HG002.m84005_220827_014912_s1.GRCh38.deepvariant.phased.vcf.gz "$dir"/HG002_ont_DUPLEX_hg38_sorted_all.bam \
		  --ignore-read-groups 

#!/bin/bash -l

#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J straglr
#SBATCH -t 100:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/straglr/log/gt-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/straglr/log/gt-%j.err-%N
#SBATCH -N 1
#SBATCH -V
#SBATCH -n 1
#SBATCH -c 3
#SBATCH --export ALL

source ~/.bashrc
conda deactivate
module purge
module load anaconda3/2021.05-xvr4x2l
conda activate straglr
which python

cd /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/straglr
ref="/tscc/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"
chr="$1"

python ./straglr/straglr.py /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG002.m84005_220827_014912_s1.GRCh38.bam $ref /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/straglr/output/pacbio_HG002_chr"$chr"_apr23 \
		--loci /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/metadata/GRCh38.hipstr_reference_pre_chr"$chr"_straglr.bed \
		--chroms chr"$chr" \
	       	--min_support 4 \
		--tmpdir /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/straglr/tmp \
	       	--genotype_in_size \
	        --nprocs 1 \
		--min_str_len 1 \
		--max_str_len 1000




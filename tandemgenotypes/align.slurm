#!/bin/bash -l

#SBATCH -p condo
#SBATCH -q condo
#SBATCH --account=ddp268
#SBATCH -J tandemgenotypes
#SBATCH -t 100:00:00
#SBATCH -o /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/tandemgenotypes/log/gt-%j.err-%N
#SBATCH -e /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/tandemgenotypes/log/gt-%j.err-%N
#SBATCH -N 1
#SBATCH -V
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --export ALL

module load samtools

dir="/tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/tandemgenotypes"
ref="/tscc/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"

time samtools view -bh /tscc/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG002.m84005_220827_014912_s1.GRCh38.bam chr21 | samtools fastq > $dir/HG002_chr21_HiFi.fastq &&

time /tscc/nfs/home/hziaeija/bin/lastdb -P20 -uRY4 "$dir"/mydb $ref &&

time /tscc/nfs/home/hziaeija/bin/last-train -P20 -Q0 "$dir"/mydb $dir/HG002_chr21_HiFi.fastq  > "$dir"/HG002_chr21_HiFi.par &&

time /tscc/nfs/home/hziaeija/bin/lastal -P20 --split -p "$dir"/HG002_chr21_HiFi.par \
"$dir"/mydb \
$dir/HG002_chr21_HiFi.fastq | gzip > "$dir"/HG002_chr21_HiFi.maf.gz

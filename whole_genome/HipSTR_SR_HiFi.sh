#PBS -q hotel 
#PBS -N HipSTR_SR_HiFi
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


dir=/projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome

time /projects/ps-gymreklab/helia/GIAB_meeting/strawman_benchmark/hipstr/HipSTR/HipSTR \
--bams /projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads/HG002.m84005_220827_014912_s1.GRCh38.bam \
--fasta /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
--regions "$dir"/metadata/GRCh38.hipstr_reference_pre.bed \
--str-vcf "$dir"/output/HG002_wg_HiFi_HipSTR_short_chr"$chr".vcf.gz \
--bam-samps HG002 --bam-libs HG002 \
--def-stutter-model \
--max-str-len 1000 \
--chrom chr"$chr" \
--max-flank-indel 1 \
--use-unpaired \
--no-rmdup \
--min-reads 5 \
--haploid-chrs chrX,chrY \
--output-filters \
--log "$dir"/short_HiFi_chr"$chr".txt

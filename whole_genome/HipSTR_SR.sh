#PBS -q hotel 
#PBS -N HipSTR_SR
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


dir=/projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome

time /projects/ps-gymreklab/helia/GIAB_meeting/strawman_benchmark/hipstr/HipSTR/HipSTR \
--bams "$dir"/Illumina_reads/HG002.GRCh38.2x250.bam \
--fasta /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
--regions "$dir"/metadata/GRCh38.hipstr_reference_pre.bed \
--str-vcf "$dir"/output/HG002_wg_Illumina_HipSTR_short.vcf.gz \
--bam-samps HG002 --bam-libs HG002 \
--max-str-len 10000 \
--def-stutter-model \
--max-flank-indel 1 \
--min-reads 5 \
--haploid-chrs chrX,chrY \
--output-filters

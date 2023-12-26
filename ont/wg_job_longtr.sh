#PBS -q hotel 
#PBS -N longtr_hg002
#PBS -l nodes=1:ppn=10:cascade
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group

cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"

time /projects/ps-gymreklab/helia/HipSTR_LR/LongSTR/LongTR \
--bams HG002_ont_DUPLEX_hg38_sorted_all.bam \
--fasta $ref \
--regions /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/hipstr_adotto_set/hipstr_catalog_adotto.bed \
--tr-vcf output/HG002_wg_longtr_ont_dec14_chr"$chr".vcf.gz \
--bam-samps HG002 --bam-libs HG002 \
--min-reads 4 \
--chrom chr"$chr" \
--output-filters \
--haploid-chrs chrX,chrY \
--max-tr-len 10000 \
--min-sum-qual -1e18 \
--skip-assembly \
--phased-bam \
--indel-flank-len 25 \
--log log/longtr_dec14_chr"$chr".txt

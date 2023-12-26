#PBS -q hotel 
#PBS -N trio_hipstr
#PBS -l nodes=1:ppn=1:cascade
#PBS -l walltime=20:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/homopolymer_fix/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/homopolymer_fix/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/homopolymer_fix
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"
read_addr="/projects/ps-gymreklab/helia/HipSTR_LR/tests/trio/reads"

echo $chr
/projects/ps-gymreklab/helia/HipSTR_LR/LongSTR/LongTR \
--bams "$read_addr"/HG002.m84005_220827_014912_s1.GRCh38.bam,"$read_addr"/HG003.m84010_220919_235306_s2.GRCh38.bam,"$read_addr"/HG004.m84010_220919_232145_s1.GRCh38.bam \
--fasta $ref \
--regions homopolymers.csv \
--tr-vcf LongTR_calls/trio_homo_longtr_stutter_chr"$chr".vcf.gz \
--bam-samps HG002,HG003,HG004 --bam-libs HG002,HG003,HG004 \
--min-reads 5 \
--output-filters \
--chrom chr"$chr" \
--max-tr-len 10000 \
--log log/log_longtr_homo_stutter_chr"$chr".txt \
--min-sum-qual -1e18 \
--skip-assembly \
--phased-bam \
--max-flank-indel 1 \
--stutter-in homopolymers_stutter.bed --stutter-align-len 1

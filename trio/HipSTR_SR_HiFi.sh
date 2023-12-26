#PBS -q hotel 
#PBS -N HipSTR_SR_HiFi
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/whole_genome/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group


cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/trio

time /projects/ps-gymreklab/helia/GIAB_meeting/strawman_benchmark/hipstr/HipSTR/HipSTR \
--bams reads/HG002.m84005_220827_014912_s1.GRCh38.bam,reads/HG003.m84010_220919_235306_s2.GRCh38.bam,reads/HG004.m84010_220919_232145_s1.GRCh38.bam \
--fasta /projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa \
--regions debug.bed \
--str-vcf test.vcf.gz \
--bam-samps HG002,HG003,HG004 --bam-libs HG002,HG003,HG004 \
--def-stutter-model \
--max-str-len 1000 \
--max-flank-indel 1 \
--use-unpaired \
--no-rmdup \
--10x-bams \
--min-reads 5 \
--output-filters

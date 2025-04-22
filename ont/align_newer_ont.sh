#PBS -q hotel 
#PBS -N align_new
#PBS -l nodes=1:ppn=20
#PBS -l walltime=50:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group

cd /projects/ps-gymreklab/helia/HipSTR_LR/tests/ont
ref="/projects/ps-gymreklab/helia/ensembl/hg38.analysisSet.fa"
source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

#aws s3 cp --no-sign-request s3://ont-open-data/gm24385_2023.12/hg002v1.0.fasta.gz hg002v1.0.fasta.gz
 
/projects/ps-gymreklab/helia/packages/minimap2/minimap2 -ax map-ont -t 15 $ref hg002v1.0.fasta.gz > hg002v1.0.alignd_hg38.sam &&

samtools view -bh hg002v1.0.alignd_hg38.sam -@ 20 > hg002v1.0.alignd_hg38.bam &&

samtools sort hg002v1.0.alignd_hg38.bam -@ 20 -o hg002v1.0.alignd_hg38_sorted.bam &&

samtools index hg002v1.0.alignd_hg38_sorted.bam -@ 20

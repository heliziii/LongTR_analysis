#PBS -q hotel 
#PBS -N assmebly
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/log
#PBS -e /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/log
#PBS -V
#PBS -M hziaeija@ucsd.edu
#PBS -A gymreklab-group

cd /projects/ps-gymreklab/helia/HipSTR_LR/TRGT/sequence_compare

source /projects/ps-gymreklab/helia/ensembl/venv_3.9/bin/activate

python3 assembly_comparison.py $chr ../hipstr_adotto_set/hipstr_catalog_adotto.bed > assembly_alleles_chr"$chr".csv

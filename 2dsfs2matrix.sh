#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-pca-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-pca-stderr.txt
#SBATCH -J 2poppca
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# Array variables
#x=$SLURM_ARRAY_TASK_ID

# Load zlib 1.2.8
module load zlib

#Declare other variables
sfs="PE"				# prefix of 2dsfs file to convert

##### CALCULATING GENOTYPES, GENOTYPE LIKELIHOODS, SFS, MAFS, ETC. WHILE ACCOUNTING FOR INBREEDING

date
echo "Prepping 2D SFS matrix..."
awk '{ OFS="\t"; for (i=1; i<=NF; i=i+65) {out=$i; for (j=i+1; j<i+65; j++) out=out"\t"$j; print out;};}' "$sfs".ml > "$sfs"_matrix.ml
date

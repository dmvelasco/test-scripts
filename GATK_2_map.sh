#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-GATK2map-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-GATK2map-stderr.txt
#SBATCH -J map
#SBATCH -p med
#SBATCH -a 1-67%10
#SBATCH -t 20-00:00:00
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=24000
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

##### PREREQUISITE STEPS:
# declare variables and arrays
# load zlib
# make mapping directory
# copy reference

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

##### CREATE SAMPLE PREFIX #####
# path to sample list
sample="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

# mapfile to extract sample ID and read name information, each line is array item
mapfile -s "$i" -n 1 -t id < "${sample}"
# -s number of rows to skip
# -n number of rows to read
# -t (cannot remember, but does not correspond to last item; remove leading/trailing whitespace?)
# id the array name

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
acc="${arr[0]}"

##### Declare directories
dir1="/home/dmvelasc/bin"				# program directory
dir2="/home/dmvelasc/Data/references/bwa-peach-scf"	# reference directory
dir3="/group/jrigrp3/Velasco/Prunus/fastq"		# sequence directory prefix
dir4="/group/jrigrp3/Velasco/Prunus/BAM"		# output directory
dir5="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
scratch="/scratch/dmvelasc"				# scratch directory for temporary storage

##### Load zlib 1.2.8
module load zlib

echo "create temporary directory to process sample";
date

# Make directory to copy reference
mkdir -p "$scratch"/"$acc"_temp

echo "copy reference files to temporary directory";
date

# Index reference to temp directories
cp "$dir2"/Prunus_persica_v1.0_scaffolds.fa.* "$scratch"/"$acc"_temp/
cp "$dir5"/Prunus_persica_v1.0_scaffolds.fa "$scratch"/"$acc"_temp/

echo "Map ${acc} to reference with BWA MEM";
date

# Map sample to reference
# Use below with fastq that does not need the sequence ID trimmed, i.e. non-SRR samples
srun "$dir1"/bwa mem -M -t 10 -k 10 "$scratch"/"$acc"_temp/Prunus_persica_v1.0_scaffolds.fa "$dir3"/"$acc"_1_filt.fq.gz "$dir3"/"$acc"_2_filt.fq.gz | "$dir1"/samtools view -T "$scratch"/"$acc"_temp/Prunus_persica_v1.0_scaffolds.fa - -o "$scratch"/"$acc".bam
mv "$scratch"/"$acc".bam "$dir4"/
# -t	threads
# -M	Mark shorter split hits as secondary (for Picard compatibility)
# -k	Minimum seed length [default 19]
# -r	Trigger re-seeding for a MEM longer than minSeedLen*FLOAT [1.5]; previously used 2.85 but may lead to lower accuracy

echo "remove ${scratch}/${acc}_temp directory";
date

# remove temporary reference directory recursively
rm -r "$scratch"/"$acc"_temp/

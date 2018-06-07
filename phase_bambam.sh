#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -a 1
#SBATCH --mem=96G
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

########## PHASE BAMs & EXTRACT FASTA for each CDS/GENE ##########
# Part 1:
# 3. use samtools mpileup, bcftools call, bcftools consensus
# (initial: use samtools view and samtools fasta to extract each CDS/gene FASTA sequence,
# does not quite work as expected.)
# 4. append FASTA sequence to one file with appropriate FASTA header for each
# -------------------------------------------------
# Part 2: (separate)
# 1. create multi sequence alignment for each CDS/gene with mafft
# 2. selection analysis with BUSTED

### Load modules ###
module load zlib
module load bambam

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
bin="/home/dmvelasc/bin"				# program directory

#### sample ID file
# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

echo -e "extract sample ID from Script/sample.txt using mapfile"
date

mapfile -s "$i" -n 1 -t id < "${list}"
# -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
# id is the array name (anything in this position is the array name)

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
acc="${arr[0]}"
echo -e "$acc"

# create SCRATCH DIRECTORY for temporary file placement
mkdir -p /home/dmvelasc/Projects/Prunus/Data/fasta/"$acc"_scaffolds

#### Index BAM file
echo -e "Check if HCrealign BAM file is indexed"
date

if [ ! -f "'$acc'_HCrealign.bam.bai" ]; then
  echo -e "Indexed file does not exist, indexing."
  "$bin"/samtools index "$acc"_HCrealign.bam
else
  echo -e "HCrealign BAM file index file exists, skipping to phasing."
fi

echo -e "Extract FASTA for each scaffold with hapHunt"
date

hapHunt "$acc"_HCrealign.bam
### outputs to working directory, can this be redirected? not with simple >

# move sample file directory from scratch
mv /group/jrigrp3/Velasco/Prunus/BAM/*.fasta /home/dmvelasc/Projects/Prunus/Data/fasta/"$acc"_scaffolds/

echo "hapHunt FASTA extraction from BAM file finished"
date

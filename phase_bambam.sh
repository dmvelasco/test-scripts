#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -a 9-10%1
#SBATCH --mem=32G
#SBATCH --exclude=bigmem1
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
ref="/home/Data/references/persica-SCF"

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

echo -e "Extract FASTA for scaffold with hapHunt"
date

while read z; do
  # create an array from each line
  scaffold=(`echo "$z"`)
  # declare variable for first column, which is "scaffold_#"
  chr=scaffold[0]
  # create a new BAM file with the selected scaffold
  srun "$bin"/samtools view -h -o "$acc"_"$chr".bam "$acc"_HCrealign.bam "$chr"
  # index the new BAM file
  "$bin"/samtools index "$acc"_"$chr".bam
  # phase the selected scaffold
  hapHunt "$acc"_"$chr".bam
  ### outputs to working directory, not easily redirected
  # move and rename file
  mv "$chr".fasta /home/dmvelasc/Projects/Prunus/Data/fasta/"$acc"_"$chr".fa
  # remove selected scaffold and associated index file
  rm "$acc"_"$chr".bam "$acc"_"$chr".bam.bai
done < "$ref"/scaffolds.txt

echo "hapHunt FASTA extraction from BAM file finished"
date

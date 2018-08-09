#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fasta/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-mafft_prep.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-mafft_prep.txt
#SBATCH -J fastacat
#SBATCH -p bigmemh
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 2-00:00:00
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# Load zlib 1.2.8
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"				# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"	# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
dir4="/scratch/dmvelasc/fasta-concat"			# scratch directory
dir5="/group/jrigrp3/Velasco/Prunus/fasta"		# directory of CDS fasta sequences

# concatenate fasta sequences from each sample for each gene
# output multi-sequence fasta for use in MAFFT multi-sequence alignment program

# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

####################
### Begin script ###
####################
# create scratch directory for temporary file placement
mkdir -p /scratch/dmvelasc/fasta-concat/

# create multi-sequence FASTA for each gene and CDS by concatenating FASTAs by ID from each sample
# Three steps:
# 1. for loop establishes sample ID
# 2. if statement determines if the sample directory with the fasta files exists
# and performs concatentation steps
# 3. while loop goes through each gene

echo "begin CDS FASTA concatenation"
date

for i in {0..66}; do
  mapfile -s "$i" -n 1 -t id < "${list}"
  # -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
  # id is the array name (anything in this position is the array name)

  # create an array from each two column line
  arr=( `echo "${id[0]}"` )

  # declare variables, created from array
  acc="${arr[0]}"
  echo -e "$acc"

  if [ -d "$acc" ]; then
    while read p; do
         # concatenate the fasta files for each sample by looping through the array
           # gene FASTA concatenation
           cat "$dir5"/"$acc"/"$p"_"$acc"_gene.fa >> "$dir4"/"$p"_gene.fa
           echo -e "\n" >> "$dir4"/"$p"_gene.fa
           # cds FASTA concatenation
           cat "$dir5"/"$acc"/"$p"_"$acc"_cds.fa >> "$dir4"/"$p"_cds.fa
           echo -e "\n" >> "$dir4"/"$p"_cds.fa
    done < "$dir3"/Prunus_persica_v1.0_genes_list.gff3
  fi
done

echo "end CDS FASTA concatenation"
date

# move files and remove directory
echo "move concatenated CDS FASTA files and remove temporary fasta scratch directory"
date
mv "$dir4" "$dir5"/


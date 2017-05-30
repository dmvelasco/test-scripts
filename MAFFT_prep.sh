#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-mafft.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-mafft.txt
#SBATCH -J fastcat
#SBATCH -p med
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 2-00:00:00
set -e
set -u

# Load zlib 1.2.8
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/scratch/dmvelasc/fasta"					# scratch directory
dir5="/home/dmvelasc/Projects/Prunus/Analysis/genetree"		# directory of CDS fasta sequences
								# each sample has a separate directory

# concatenate fasta sequences from each sample for each gene
# output multi-sequence fasta for use in MAFFT multi-sequence alignment program

# full set of initial samples
declare -a id=(PS02 PD01 PP15)

####################
### Begin script ###
####################
echo "begin CDS FASTA concatenation"
date

# create scratch directory for temporary file placement
mkdir -p /scratch/dmvelasc/fasta/

# create multi-sequence FASTA for each gene CDS FASTA
# while loop goes through each gene in the list file (end of loop)
while read p; do
     touch "$dir4"/"$p".fa
     # concatenate the fasta files for each sample by looping throug the array
     for each in "${id[@]}"
     do
        cat "$dir5"/"$each"/"$p"_"$each".fa >> "$dir4"/"$p".fa
        echo -e "\n" >> "$dir4"/"$p".fa
     done
done < "$dir3"/Prunus_persica_v1.0_genes_list.gff3

echo "end CDS FASTA concatenation"
date

# move files and remove directory
echo "move concatenated CDS FASTA files and remove temporary fasta scratch directory"
date

mv "$dir4"/*.fa "$dir5"/
rm -rf "$dir4"

#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-mafft.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-mafft.txt
#SBATCH -J mafft
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2-00:00:00
set -e
set -u

# Load zlib 1.2.8
module load zlib
#PAML is not a module on farm

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"


while read p; do
        # create multi-sequence fasta for each gene
        touch "$dir4"/"$p".fa
        for i in "${id[@]}"
        do
           cat "$dir5"/"${id[$i]}"/"$p"_"${id[$i]}" >> "$dir4"/"$p".fa
        done
done < "$dir3"/Prunus_persica_v1.0_genes_list.gff3

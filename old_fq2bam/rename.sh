#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /group/jrigrp3/Velasco/Prunus/BAM2/%j-stdout-rename.txt
#SBATCH -e /group/jrigrp3/Velasco/Prunus/BAM2/%j-stderr-rename.txt
#SBATCH -J rename
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u


# %j is job allocation number

# removing _scf.sort.bam extension and replacing with .bam

# Delcare directories
dir1="/group/jrigrp3/Velasco/Prunus/BAM2"

for i in $dir1/*.bam; do
	file="${i##*/}" 			# gathers file name without path
	id="${file%_*}" 			# retain only core id
	mv "$id"_scf.sort.bam "$id".bam		# move to new filename
done

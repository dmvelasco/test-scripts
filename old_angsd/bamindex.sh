#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stats-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stats-stderr.txt
#SBATCH -J index
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Declare program directory
dir1="/home/dmvelasc/bin"

# Load zlib 1.2.8
module load zlib

# create index files (.bai) for each BAM file

for i in *.bam; do
        file="${i##*/}" #gathers the file name with extension
	"$dir1"/samtools index "$file"
done

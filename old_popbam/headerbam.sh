#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /group/jrigrp3/Velasco/Prunus/BAM2/%j-merge-stdout.txt
#SBATCH -e /group/jrigrp3/Velasco/Prunus/BAM2/%j-merge-stderr.txt
#SBATCH -J merge
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u

# Declare directories
dir1="/home/dmvelasc/bin"                               # program directory
dir2="/group/jrigrp3/Velasco/Prunus/Analysis"		# final directory

# SAMtools merge of BAM files for scaffold 1
"$dir1"/samtools reheader all.newheader.sam "$dir2"/all.prunus.bam > "$dir2"/allnew.prunus.bam
rm "$dir2"/all.prunus.bam
mv "$dir2"/allnew.prunus.bam "$dir2"/all.prunus.bam
"$dir1"/samtools index "$dir2"/all.prunus.bam

#In case need to change only the header, in the future try samtools reheader
#usage: samtools reheader <in.header.sam> <in.bam>
#"Replace the header in in.bam with the header in in.header.sam. This command is much faster than replacing the header with a BAM->SAM->BAM conversion."
#better than the 20 hours it takes to remerge all bam files

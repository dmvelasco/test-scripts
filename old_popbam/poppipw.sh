#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /group/jrigrp3/Velasco/Prunus/BAM2/%j-poppi-stdout.txt
#SBATCH -e /group/jrigrp3/Velasco/Prunus/BAM2/%j-poppi-stderr.txt
#SBATCH -J poppi
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 16
set -e
set -u


# Declare directories
dir1="/home/dmvelasc/bin"                               # program directory
dir2="/home/dmvelasc/Data/references/persica-SCF"       # FASTA reference directory
dir3="/group/jrigrp3/Velasco/Prunus/Analysis"           # final directory

# SAMtools merge of BAM files for scaffold 1
"$dir1"/samtools merge -r -h all.header.sam "$dir3"/all.prunus.bam -b all.prunusbam.txt
"$dir1"/samtools index "$dir3"/all.prunus.bam

# copy peach reference scaffolds fasta
cp "$dir2"/Prunus_persica_v1.0_scaffolds.fa "$dir3"/   # copy reference to final directory

# POPBAM pi along genome (scaffold 1)
for i in {1..8}; do
"$dir1"/popbam nucdiv -w 1 -f Prunus_persica_v1.0_scaffolds.fa -a 30 all.prunus.bam scaffold_"$i" > all.prunus-1kb_"$i".pi
done

#         -w  INT     use sliding window of size (kb)
#         -n  INT     minimum sample size per population   [ default: all samples present ]
#         -f  FILE    Reference fastA file
#         -a  INT     minimum map quality                  [ default: 13 ]

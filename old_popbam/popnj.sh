#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-popnj-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-popnj-stderr.txt
#SBATCH -J popnj
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u

# Declare directories
bin="/home/dmvelasc/bin"                               # program directory

# ORIGINAL PARAMETERS
#for i in {1..8}; do
#"$bin"/popbam tree -f Prunus_persica_v1.0_scaffolds.fa /all.prunus.bam scaffold_"$i" > all.prunus_"$i".tree
#done

for i in {1..8}; do
"$bin"/popbam tree -f Prunus_persica_v1.0_scaffolds.fa -a 20 -b 20 all.prunus.bam scaffold_"$i" > all.prunus_"$i".20.tree
done

# Usage:   popbam tree [options] <in.bam> [region]

# Options:
#	  -i          base qualities are Illumina 1.3+     [ default: Sanger ]
#         -h  FILE    Input header file                    [ default: none ]
#         -d  STR     distance (pdist or jc)               [ default: pdist ]
#         -w  INT     use sliding window of size (kb)
#         -k  INT     minimum number of sites in window    [ default: 10 ]
#         -f  FILE    Reference fastA file
#         -m  INT     minimum read coverage                [ default: 3 ]
#         -x  INT     maximum read coverage                [ default: 255 ]
#         -q  INT     minimum rms mapping quality          [ default: 25 ]
#         -s  INT     minimum snp quality                  [ default: 25 ]
#         -a  INT     minimum map quality                  [ default: 13 ]
#         -b  INT     minimum base quality                 [ default: 13 ]

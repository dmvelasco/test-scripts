#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /group/jrigrp3/Velasco/Prunus/BAM2/%j-popsfs-stdout.txt
#SBATCH -e /group/jrigrp3/Velasco/Prunus/BAM2/%j-popsfs-stderr.txt
#SBATCH -J popsfs
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u


# Declare directories
dir1="/home/dmvelasc/bin"                               # program directory

# POPBAM pi along genome (scaffold 1)
for i in {1..8}; do
"$dir1"/popbam sfs -w 1 -p PC01 -f Prunus_persica_v1.0_scaffolds.fa -a 30 all.prunus.bam scaffold_"$i" > all.prunus-1kb_"$i".sfs
done

#Usage:   popbam sfs [options] <in.bam> [region]

#Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]
#         -h  FILE    Input header file                    [ default: none ]
#         -w  INT     use sliding window of size (kb)
#         -p  STR     sample name of outgroup              [ default: reference ]
#         -f  FILE    Reference fastA file
#         -m  INT     minimum read coverage                [ default: 3 ]
#         -x  INT     maximum read coverage                [ default: 255 ]
#         -q  INT     minimum rms mapping quality          [ default: 25 ]
#         -s  INT     minimum snp quality                  [ default: 25 ]
#         -a  INT     minimum map quality                  [ default: 13 ]
#         -b  INT     minimum base quality                 [ default: 13 ]

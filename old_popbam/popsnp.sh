#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-popsnp-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-popsnp-stderr.txt
#SBATCH -J popsnp
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u

# Declare directories
dir1="/home/dmvelasc/bin"                               # program directory

# POPBAM snps by scaffold
for i in {1..8}; do
"$dir1"/popbam snp -p PC01 -o 1 -f Prunus_persica_v1.0_scaffolds.fa -a 30 -z 0.05 all.prunus.bam scaffold_"$i" > all.prunus.MQ30_scaffold"$i"_z05.snp
done

for i in {1..8}; do
"$dir1"/popbam snp -p PC01 -o 1 -f Prunus_persica_v1.0_scaffolds.fa -a 30 -z 0.5 all.prunus.bam scaffold_"$i" > all.prunus.MQ30_scaffold"$i"_z50.snp
done

#did not work with following commands included, removing temporarily 2014-12-01, add back in singly to test which is/are causing hang up
# -v <- see note below
# -z 0.05 or 0.5
# -a 30 <- restored
# -b 20
# -w 1 <- see note below

# Usage:   popbam snp [options] <in.bam> [region]

# Options:
#	  -i          base qualities are Illumina 1.3+     [ default: Sanger ]
#         -h  FILE    Input header file                    [ default: none ]
#         -v          output variant sites only            [ default: All sites ]
#         -z  FLT     output heterozygous base calls       [ default: Consensus ] FLT=float, default 0.0001
# based on Verde et al. paper the nucleotide diversity for P. davidiana is 4.8 x 10^-3
#         -w  INT     use sliding window of size (kb)
#         -p  STR     sample name of outgroup              [ default: reference ]
#         -o  INT     output format                        [ default: 0 ]
#                     0 : popbam snp format
#                     1 : SweepFinder snp format
#                     2 : MS format
#         -f  FILE    Reference fastA file
#         -m  INT     minimum read coverage                [ default: 3 ]
#         -x  INT     maximum read coverage                [ default: 255 ]
#         -q  INT     minimum rms mapping quality          [ default: 25 ]
#         -s  INT     minimum snp quality                  [ default: 25 ]
#         -a  INT     minimum map quality                  [ default: 13 ]
#         -b  INT     minimum base quality                 [ default: 13 ]

#Garrigan notes from github
#-w
#intended to be used with the ms output option
#each window is treated as a separate sample in a typical ms run
#no effect when invoked in the other two output modes

#-v
#only available when the output mode is set to popbam native snp format
#allows users to count the number of aligned sites

#-h
#used to call heterozygous sites and only works when the output is in native popbam format
#otherwise the consensus base is used

#-p
#allows users to specify the outgroup sample to use
#by default the reference is assumed to be the outgroup
#The name of the outgroup sample must match the `SM' field value in the `RG' header in the input BAM file.
#Specification of the outgroup sequence is used to polarize ancestral and derived states of polymorphic sites,
#and thus its effects can only be seen when output is in SweepFinder or ms format.

#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis/bcf
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcfconvert.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcfconvert.txt
#SBATCH -J vcf
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u


# %j is job allocation number

# load zlib 1.2.8 otherwise uses zilb 1.2.3
module load zlib

# Declare directories
dir1="/home/dmvelasc/Software" # directory with perl script for converting vcf to hapmap format
dir2="/group/jrigrp3/Velasco/Prunus/Analysis/bcf" # VCF directory prefix

cat "$dir2"/prunus.flt.vcf | "$dir1"/vcf2hp2a.pl > prunus.hmp.txt
# use zcat if gzipped file
# perl script modified because prints some extra stuff which may be incompatible with tassel
# also modified to ignore quolity... let's see if that works

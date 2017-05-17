#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf2struct.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf2struct.txt
#SBATCH -J vcf2struct
#SBATCH -p bigmeml
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Load zlib 1.2.8 not recognized
#module load zlib

# Declare directories and file prefix
dir1="/home/dmvelasc/bin"				# software binary directory
dir2="/home/dmvelasc/Software"				# general software directory
file="amygdalus"							# file prefix

##### FILTER VCF, CONVERT TO HAPMAP, RANDOMLY SELECT 10000 LOCI, CONVERT TO STRUCTURE FORMAT #####

# Filter VCF and convert to hapmap format file
"$dir1"/bcftools view -i "INFO/MQ>=30" "$file".fltID.vcf | "$dir2"/vcf2hp2.pl > "$file".hmp

# Randomly select 10K SNP loci
grep -v "^#" "$file".hmp | sort -R - | tail -10000 > "$file"10K_nohdr.hmp

# Restore header to file
grep "^rs#" "$file".hmp | cat - "$file"10K_nohdr.hmp > "$file"10K.hmp

# Convert to STRUCTURE format
# use SNPutilities interactively

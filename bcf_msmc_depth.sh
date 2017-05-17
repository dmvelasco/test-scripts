#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf.txt
#SBATCH -J vcf
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 1
set -e
set -u


# Almond and peach genotype per chromosome depth stats for Josh

# Load zlib 1.2.8
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"


#"$dir1"/bcftools index dulcis_msmc.flt.vcf.bzip
#"$dir1"/bcftools index persica_msmc.flt.vcf.bzip

touch msmc_15xdepth_perchr.txt

##### DEPTH CHECK PER CHROMOSOME #####
for i in {1..8}; do
	echo "almond scaffold "$i"\n" >> msmc_15xdepth_perchr.txt
	"$dir1"/bcftools gtcheck "$dir2"/dulcis_msmc.flt.vcf.bzip -r scaffold_"$i" >> msmc_15xdepth_perchr.txt
	echo "peach scaffold "$i"\n" >> msmc_15xdepth_perchr.txt
	"$dir1"/bcftools gtcheck "$dir2"/persica_msmc.flt.vcf.bzip -r scaffold_"$i" >> msmc_15xdepth_perchr.txt
done

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


# To compare with MSMC and Josh's SAMtools results using reduced sampling of almond and peach genotypes

# Load zlib 1.2.8
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"

##### POPGEN STATISTICS WITH VCFTOOLS #####
#"$dir1"/vcftools --gzvcf "$dir2"/dulcis_msmc20x.flt.vcf.bzip --out dulcis_msmc20x --site-pi
#"$dir1"/vcftools --gzvcf "$dir2"/dulcis_msmc20x.flt.vcf.bzip --out dulcis_msmc20x-step --window-pi 1000 --window-pi-step 50
#"$dir1"/vcftools --gzvcf "$dir2"/dulcis_msmc20x.flt.vcf.bzip --out dulcis_msmc20x --TajimaD 1000
"$dir1"/vcftools --gzvcf "$dir2"/dulcis_msmc.flt.vcf.bzip --out dulcis_msmc-step --window-pi 1000 --window-pi-step 50

#"$dir1"/vcftools --gzvcf "$dir2"/persica_msmc20x.flt.vcf.bzip --out persica_msmc20x --site-pi
#"$dir1"/vcftools --gzvcf "$dir2"/persica_msmc20x.flt.vcf.bzip --out persica_msmc20x-step --window-pi 1000 --window-pi-step 50
#"$dir1"/vcftools --gzvcf "$dir2"/persica_msmc20x.flt.vcf.bzip --out persica_msmc20x --TajimaD 1000
"$dir1"/vcftools --gzvcf "$dir2"/persica_msmc.flt.vcf.bzip --out persica_msmc-step --window-pi 1000 --window-pi-step 50

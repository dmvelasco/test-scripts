#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
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

##### FILTER SAMPLES, CALL A FEW STATISTICS #####
# pileup and filtering of genomic SNPs
# with -A removed (consider removing -v, which only calls variant sites)
#"$dir1"/samtools mpileup -BRug -t DP,SP -f "$dir3"/Prunus_persica_v1.0_scaffolds.fa -b all.prunusbam.txt | "$dir1"/bcftools call -O z -vm -f GQ,GP - > "$dir2"/prunus-A.raw.vcf.bzip
#"$dir1"/bcftools index "$dir2"/prunus-A.raw.vcf.bzip
"$dir1"/bcftools view -O z -V indels -o "$dir2"/dulcis_msmc20x.flt.vcf.bzip -s PD03.bam,PD04.bam,PD05.bam,PD06.bam,PD07.bam,PD09.bam -i "%QUAL>=30 && %AVG(DP)>20 && %MIN(GQ)>=30" "$dir2"/prunus-A.raw.vcf.bzip
"$dir1"/bcftools view -O z -V indels -o "$dir2"/persica_msmc20x.flt.vcf.bzip -s PP02.bam,PP03.bam,PP04.bam,PP05.bam,PP13.bam -i "%QUAL>=30 && %AVG(DP)>20 && %MIN(GQ)>=30" "$dir2"/prunus-A.raw.vcf.bzip

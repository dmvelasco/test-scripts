#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-vcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-vcf.txt
#SBATCH -J vcf
#SBATCH -p bigmemm
#SBATCH -a 1-6
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u


# %A is array job ID
# %a is array job index
# %j is job allocation number

x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Load zlib 1.2.8
module load zlib

# Declare BAM file list array
declare -a group=(dulcis persica almond peach amygdalus prunus)

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Script/old_samtools"	# script directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory

##### PILEUP, CALL VARIANTS, FILTER VARIANTS #####
# PILEUP OF SORTED, UNCOMPRESSED BAM TO BCF
# filter VCF (approximate average depth is 30X)

# pileup and filtering of genomic SNPs
"$dir1"/samtools mpileup -BRug -t DP,SP -f "$dir3"/Prunus_persica_v1.0_scaffolds.fa -b "$dir2"/all."${group["$i"]}"bam.txt | "$dir1"/bcftools call -O z -Avm -f GQ,GP - > "${group["$i"]}".raw.vcf.bzip
"$dir1"/bcftools index "${group["$i"]}".raw.vcf.bzip
"$dir1"/bcftools view -O z -o "${group["$i"]}".flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %AVG(DP)<=100 && %MIN(GQ)>=30" -r scaffold_5:10900000-12700000 "${group["$i"]}".raw.vcf.bzip

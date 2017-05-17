#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-pileup.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-pileup.txt
#SBATCH -J pileup
#SBATCH -p hi
#SBATCH -a 1-6
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u


# %A is array job ID
# %a is array job index
# %j is job allocation number

x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Loac zlib 1.2.8
module load zlib

# Declare BAM file list array
declare -a group=(dulcis persica almond peach amygdalus prunus)

# Declare directories
dir1="/home/dmvelasc/bin"				# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Script"		# output directory
dir3="/home/dmvelasc//Data/references/persica-SCF"	# FASTA reference directory

##### PILEUP, CALL VARIANTS, FILTER VARIANTS #####
# PILEUP OF SORTED, UNCOMPRESSED BAM TO BCF
# filter VCF (approximate average depth is 30X)

# peach type
"$dir1"/samtools mpileup -DBSRu -f "$dir3"/Prunus_persica_v1.0_scaffolds.fa -b "$dir2"/all."${group["$i"]}"bam.txt | "$dir1"/bcftools call -O b -Avm - > "${group["$i"]}".raw.bcf
"$dir1"/vcftools --bcf "${group["$i"]}".raw.bcf --remove-indels --minQ 30 --min-meanDP 3 --max-meanDP 100 --thin 11 --recode --filtered-sites --out "${group["$i"]}".final

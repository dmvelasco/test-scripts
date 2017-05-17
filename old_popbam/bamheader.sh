#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-bamhdr-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-bamhdr-stderr.txt
#SBATCH -J bamhdr
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Declare directories
dir1="/home/dmvelasc/bin"                               # program directory
dir2="/group/jrigrp3/Velasco/Prunus/BAM2"		# BAM directory
dir3="/group/jrigrp3/Velasco/Prunus/BAM2/test"		# output directory
# Delcare directories

for i in $dir2/PD*.bam; do
        file="${i##*/}"		# gathers file name without path
        id="${file%_*}"         # retain only core id
#	pre="${id%[0-9]*}"	# want only prefix, i.e., PD, PP, etc. <- NOT COMPLETED
#	if [ "$id" ?? "PD"  ]	# id begins with PD
#	then
		# SAMtools reheader P. dulcies BAM files; contains modified scaffold 1 length
		"$dir1"/samtools reheader angsd.header.sam "$id" > "$dir3"/"$id"
		"$dir1"/samtools index "$dir3"/"$id"
#	fi
done

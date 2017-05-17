#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2/test
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-map-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-map-stderr.txt
#SBATCH -J bwa2bam
#SBATCH -p serial
#SBATCH -a 1-8
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u

# %A is array job ID
# %a is array job index

##### PREREQUISITE STEPS: declare variables and arrays (check), load zlib (check), make mapping directory (check), copy reference (check)

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare arrays (change depending on data organization)
declare -a id=(PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
prefix="${id["$i"]}"

# make directory to copy reference
mkdir "$prefix"_temp

# Declare directories
dir1="/home/dmvelasc/bin"				# program directory
dir2="/home/dmvelasc/Data/references/bwa-peach2"	# reference directory
dir3="/group/jrigrp3/Velasco/Prunus/sickle"		# sequence directory prefix
dir4="/group/jrigrp3/Velasco/Prunus/BAM2/test"		# output directory
dir5="/home/dmvelasc/Data/references/bwa-peach2"	# FASTA reference directory

# Declare other variables
genome="Prunus_persica_v2.0.a1_scaffolds.fasta"

# Load zlib 1.2.8
module load zlib

# Index reference to temp directories
cp "$dir2"/"$genome"* "$prefix"_temp/

##### Map to reference with BWA MEM #####

"$dir1"/bwa mem -t 8 -k 10 -r 2.85 "$prefix"_temp/"$genome" "$dir3"/sickle_"$prefix"_1.fq.gz "$dir3"/sickle_"$prefix"_2.fq.gz | "$dir1"/samtools view -uT "$prefix"_temp/"$genome" - -o "$dir4"/"$prefix".bam

# remove temporary reference directory recursively
rm -r "$prefix"_temp/


##### SAM AND BAM CONVERSIONS AND SORTING #####

# sort BAM file - CONFIRM MEMORY QTY
# -@ is number of sorting and compression threads; -m is memory per thread
# new options v 1.3
# -T temporary file prefix, -o final output file, -l compression level (0=uncompressed, 9= best)
#"$dir1"/samtools sort -@ 8 -m 1G -l 0 -T "$dir4"/"$prefix".sorted -O BAM -o "$dir4"/"$prefix".sort.bam "$dir4"/"$prefix".bam
"$dir1"/samtools sort -@ 8 -m 1G -l 0 -T "$prefix".sorted -O BAM -o "$prefix".sort.bam "$prefix".bam

# remove unsorted BAM
rm "$prefix".bam

# convert back to uncompressed BAM
#"$dir1"/samtools view -bhu "$dir4"/"$prefix".sorted.bam -o "$dir4"/"$prefix".sort.bam

# rename sorted BAM to <prefix>.bam
mv "$prefix".sort.bam "$prefix".bam

# samtools rmdup to remove duplicate reads should work with v 1.3, test later

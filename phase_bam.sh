#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-map-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-map-stderr.txt
#SBATCH -J bwa2bam
#SBATCH -p serial
#SBATCH -a 1-20%5
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=24000
set -e
set -u


########## PHASE AND EXTRACT CDS AND/OR GENE FASTA FILES ##########
# Step 1: prep prefix information to cycle through each BAM file
# (using the GATK produced HCrealign.bam files)
# Step 2: use phase to output at least two new phased BAM files
# Step 3: use view and fasta to extract each CDS/gene FASTA sequence
# Step 4: append FASTA sequences to one file with  appropriate FASTA header for each
#---------------------------------------
# Step 5: create multi sequence alignment for each CDS/gene with mafft
# Step 6: selection analysis

##### PREREQUISITE STEPS #####
# declare variables and arrays
# load zlib
# make mapping directory
# copy reference

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare prefix array
declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)

acc="${abbr["$i"]}"

# make directory to copy reference
#mkdir "$acc"_temp

# Declare directories
dir1="/home/dmvelasc/bin"				# program directory
BAMdir="/group/jrigrp3/Velasco/Prunus/BAM"		# BAM directory
dir2="/home/dmvelasc/Data/references/bwa-peach-scf"	# reference directory
dir3="/home/dmvelasc/Projects/Prunus/Data/fastq"	# sequence directory prefix
dir4="/home/dmvelasc/Projects/Prunus/Data/BAM"		# output directory
dir5="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory

# Load zlib 1.2.8
module load zlib

# Index reference to temp directories
cp "$dir2"/Prunus_persica_v1.0_scaffolds.fa.* "$acc"_temp/
cp "$dir5"/Prunus_persica_v1.0_scaffolds.fa "$acc"_temp/

##### Map to reference with BWA MEM #####

srun "$dir1"/bwa mem -M -t 14 -k 10 "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa "$dir3"/"$acc"_1_filt.fq.gz "$dir3"/"$acc"_2_filt.fq.gz | "$dir1"/samtools view -T "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa - -o "$dir4"/"$acc".bam
# -t	threads
# -M	Mark shorter split hits as secondary (for Picard compatibility)
# -k	Minimum seed length [default 19]
# -r	Trigger re-seeding for a MEM longer than minSeedLen*FLOAT [1.5] have previously used 2.85 but may lead to lower accuracy

# remove temporary reference directory recursively
rm -r "$acc"_temp/


##### SAM AND BAM CONVERSIONS, SORTING, AND CLEAN UP #####

# sort BAM file - CONFIRM MEMORY QTY
# -@ is number of sorting and compression threads; -m is memory per thread
srun "$dir1"/samtools sort -l 9 -@ 14 -m 1G -o "$dir4"/"$acc".sorted.bam -T "$dir4"/"$acc".sorting "$dir4"/"$acc".bam

# remove unsorted BAM
rm "$dir4"/"$acc".bam

# convert back to uncompressed BAM
"$dir1"/samtools view -bh "$dir4"/"$acc".sorted.bam -o "$dir4"/"$acc".bam

# remove sorted compressed BAM, final BAM file is <prefix>.bam
rm "$dir4"/"$acc".sorted.bam


##### REMOVE DUPLICATE READS #####

# remove duplicate reads
"$dir1"/samtools rmdup -S "$dir4"/"$acc".bam "$dir4"/"$acc".nodup.bam

# remove BAM with duplicate reads
rm "$dir4"/"$acc".bam

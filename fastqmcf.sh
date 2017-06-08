#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fastq/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-fastqmcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-fastqmcf.txt
#SBATCH -a 1-67%10
#SBATCH -J fastqmcf
#SBATCH -p med
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 10:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

############################################################################
### fastqmcf is a fast QC program that trims adapter sequences           ###
### and low quality bases from FASTQ files                               ###
### adapter trimming utilizes a supplied FASTA file of adapter sequences ###
############################################################################

# load modules
module load fastqmcf

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare variables
# path to FASTA adapter file
adapter="/home/dmvelasc/Projects/Prunus/Script/adapter.fa"
# path to sample list
sample="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

# mapfile to extract sample ID and read name information, each line is array item
mapfile -s "$i" -n 1 -t id < "${sample}"
# -s number of rows to skip
# -n number of rows to read
# -t (cannot remember, but does not correspond to last item)
# id the array name

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
reads="${arr[1]}"
acc="${arr[0]}"

# run fastq-mcf to trim reads
fastq-mcf -t16 "$adapter" "$reads"_1.fastq.gz -o "$acc"_1_filt.fq.gz "$reads"_2.fastq.gz  -o "$acc"_2_filt.fq.gz

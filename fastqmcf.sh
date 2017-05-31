#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fastq/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-fastqmcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-fastqmcf.txt
#SBATCH -a 1-67
#SBATCH -J fastqmcf
#SBATCH -p bigmemm
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
adapter="/home/dmvelasc/Projects/Prunus/Script/adapter.fa
# path to sample list
sample="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

# Declare prefix array <- if mapfile below works for array jobs
#declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210 SRR068360 SRR068361 SRR502982 SRR502983 SRR502984 SRR502985 SRR502986 SRR502987 SRR502990 SRR502992 SRR502994 SRR502995 SRR502997 SRR502998 SRR503000 SRR765679 SRR765838 SRR765850 SRR765861 SRR3237762 SRR3138171 SRR3141016 SRR3141018 SRR3138115 SRR3138121 SRR3138123 SRR3138169 SRR3237746 SRR3141019 SRR3136174 SRR3136179 SRR3136181 SRR3136183 SRR765861 SRR3138129 SRR3138145 SRR3138147 SRR3141049 SRR3141073 SRR3141113 SRR3141248 SRR3138117 SRR3138168 SRR3138132 SRR3141238 SRR3141181)
#declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10 PP14 PP12 PV01 PP08 PS01 PP01 PP06 PP05 PP04 PP03 PP02 PP10 PP09 PG01 PP07 PD14 PD13 PD12 PD11 PV03 PV04 PV05 PV06 PG02 PG04 PG05 PS04 PM01 PM02 PM03 PM04 PM05 PM06 PD11 PP37 PP39 PP40 PD16 PD17 PD18 PD21 PG03 PS03 PP38 PD20 PD19)

#reads="${main["$i"]}"
#acc="${abbr["$i"]}"

mapfile -s "$i" -n 1 -t id < "${sample}"
arr=(`echo ${id[0]}`)
reads="${arr[0]}"
acc="${arr[1]}"

fastq-mcf -t16 "$adapter" "$reads"_1.fastq.gz -o "$acc"_1_filt.fq.gz "$reads"_2.fastq.gz  -o "$acc"_2_filt.fq.gz

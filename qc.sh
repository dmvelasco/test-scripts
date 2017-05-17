#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Peach_GDR
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/qc-stdout-%A_%a.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/qc-stderr-%A_%a.txt
#SBATCH -J qc
#SBATCH -p serial
#SBATCH -a 1-23
#SBATCH -t 3-00:00:00
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u


module load sickle
module load zlib

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare prefix array
declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210)
declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
#declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502985 SRR502994 SRR502992 SRR502990 SRR502987 SRR502986 SRR503000 SRR502983 SRR502997 SRR502983 SRR502997 SRR502995 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
#declare -a pub=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES
declare -a adapt=()

reads="${main["$i"]}"
acc="${abbr["$i"]}"

# Declare directories (also change depending on organization)
dir1="/home/dmvelasc/Software/scythe" # program directory
dir2="/home/dmvelasc/Projects/Prunus/Data/fastq" # sequence directory prefix
dir3="/home/dmvelasc/Projects/Prunus/Analysis" # outfile directory prefix


## scythe for adapter trimming (ensure illumina_adapters are correct)
"$dir1"/scythe -a "$dir1"/illumina_adapters.fa -o "$dir2"/scythe_"$acc"_1.fq.gz "$dir2"/"$reads"_1.fq.gz
"$dir1"/scythe -a "$dir1"/illumina_adapters.fa -o "$dir2"/scythe_"$acc"_2.fq.gz "$dir2"/"$reads"_2.fq.gz

## sickle for low quality trimming
sickle pe -f "$dir2"/scythe_"$acc"_1.fq.gz -r "$dir2"/scythe_"$acc"_2.fq.gz -t sanger -o "$dir2"/sickle_"$acc"_1.fq -p "$dir2"/sickle_"$acc"_2.fq -s "$dir2"/sickle_"$acc"_singles.fq

gzip "$dir2"/sickle_"$acc"_1.fq
gzip "$dir2"/sickle_"$acc"_2.fq
gzip "$dir2"/sickle_"$acc"_singles.fq

rm "$dir2"/scythe_"$acc"_1.fq.gz "$dir2"/scythe_"$acc"_2.fq.gz

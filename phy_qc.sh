#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-qc-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-qc-stderr.txt
#SBATCH -a 1-8
#SBATCH -J qc
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u


module load sickle
module load zlib

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$x-1

## Declare arrays (these will need to change depending on sequence data organization)
## most sequenced genotypes (no special adapter sequences for scythe)
#declare -a id=(PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)

# UCB sequenced almonds
declare -a id=(PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
declare -a truseq=(CGATGT TGACCA ACAGTG GCCAAT CAGATC CTTGTA AGTCAA AGTTCC)

# Declare directories (also change depending on organization)
dir1="/home/dmvelasc/Software/scythe" # program directory
dir2="/group/jrigrp3/Velasco/Prunus/fastq" # sequence directory prefix
dir3="/group/jrigrp3/Velasco/Prunus/sickle" # outfile directory prefix


### SCYTHE
## adapter trimming

## standard illumina adapters
#"$dir1"/scythe -a "$dir1"/illumina_adapters.fa -o "$dir2"/scythe_"${id["$i"]}"_1.fq.gz "$dir2"/"${id["$i"]}"_1.fq.gz
#"$dir1"/scythe -a "$dir1"/illumina_adapters.fa -o "$dir2"/scythe_"${id["$i"]}"_2.fq.gz "$dir2"/"${id["$i"]}"_2.fq.gz

## truseq adapters (specific to UCB sequenced almonds)
"$dir1"/scythe -a "$dir1"/truseq_adapters_"${truseq["$i"]}".fasta -o "$dir2"/scythe_"${id["$i"]}"_1.fq.gz "$dir2"/"${id["$i"]}"_R1.gz
"$dir1"/scythe -a "$dir1"/truseq_adapters_"${truseq["$i"]}".fasta -o "$dir2"/scythe_"${id["$i"]}"_2.fq.gz "$dir2"/"${id["$i"]}"_R2.gz

## sickle for low quality trimming
sickle pe -f "$dir2"/scythe_"${id["$i"]}"_1.fq.gz -r "$dir2"/scythe_"${id["$i"]}"_2.fq.gz -t sanger -o "$dir2"/sickle_"${id["$i"]}"_1.fq -p "$dir2"/sickle_"${id["$i"]}"_2.fq -s "$dir2"/sickle_"${id["$i"]}"_singles.fq

gzip -c "$dir2"/sickle_"${id["$i"]}"_1.fq > "$dir2"/sickle_"${id["$i"]}"_1.fq.gz
gzip -c "$dir2"/sickle_"${id["$i"]}"_2.fq > "$dir2"/sickle_"${id["$i"]}"_2.fq.gz
gzip -c "$dir2"/sickle_"${id["$i"]}"_singles.fq > "$dir2"/sickle_"${id["$i"]}"_singles.fq.gz
rm "$dir2"/sickle_"${id["$i"]}"_1.fq "$dir2"/sickle_"${id["$i"]}"_2.fq "$dir2"/sickle_"${id["$i"]}"_singles.fq

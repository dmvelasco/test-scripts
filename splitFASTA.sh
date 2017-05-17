#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/fastq
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-fastqmcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-fastqmcf.txt
#SBATCH -a 1-20%5
#SBATCH -J fastqmcf
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 6:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
prefix="${id["$i"]}"

#Simple UNIX code for splitting FASTA on '>'
#see https://www.biostars.org/p/2226/
#Solutions offered by users Biomonika (Noolean) and new

csplit -f "$prefix"fa_ -s file.fa '/>/' {*}

#for testing
#for a in prefix_*; do echo $a; mv $a $(head -1 $a | tr -d '>' | tr " " "_" | tr ":" "_")_prefix.fa; done;

#for production (eliminate echo statement)
for a in prefix_*; do mv $a $(head -1 $a | tr -d '>' | tr " " "_" | tr ":" "_")_prefix.fa; done;

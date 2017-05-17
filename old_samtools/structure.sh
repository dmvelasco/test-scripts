#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-structure.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-structure.txt
#SBATCH -J structure
#SBATCH -p serial
#SBATCH -a 1-6
#SBATCH -n 1
#SBATCH -c 4

set -e
set -u

# %A is array job ID
# %a is array job index

w=$SLURM_ARRAY_TASK_ID
x=$(( w+6 )) # where w = $SLURM_ARRAY_TASK_ID

# Declare directories
dir1="/home/dmvelasc/bin"	# program binary location
dir2="/home/dmvelasc/Projects/Prunus/Analysis/STRUCTURE/2014-11-19/results"

for i in {1..20}; do
"$dir1"/structure -K "$x" -o "$dir2"/Prunus_K"$x"_"$i"
done


#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree/fasta
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-mafft_prep.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-mafft_prep.txt
#SBATCH -J mafft
#SBATCH -a 1-27585%10
#SBATCH -p bigmemh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

# eventual number of cycles
# -a 1-27585%50

module load mafft

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

declare -a id=(PS02 PD01 PP15)
#declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# need to modify so the array is pulled from a file with list of IDs
# also make rest of this script into a bash script that is just called by this slurm script

# NEED array for the input

##### STEP 1 #####
# align sequences
linsi input > output

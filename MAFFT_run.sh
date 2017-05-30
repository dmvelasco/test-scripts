#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree/fasta
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-mafft_prep.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-mafft_prep.txt
#SBATCH -J mafft
#SBATCH -a 1-10%2
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 20:00:00
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

# Declare directories
dir1="/home/dmvelasc/bin"                                       # software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/genetree"              # VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"               # FASTA reference directory

# basic set up
# input is multi-sequence fasta

##### STEP 1: CREATE ARRARY OF GENE IDs #####
# creates array from gene ID file, each line (gene) a separate array item
# can use to create array job
# 27864 genes in gene ID list file

mapfile -t gene < "$dir3"/Prunus_persica_v1.0_genes_list.gff3


##### STEP 2: ALIGN MULTI FASTA SEQUENCE #####
mafft --localpair --maxiterate 1000 --phylipout "$dir2"/"${gene["$i"]}".fa > "$dir2"/"${gene["$i"]}"_aln.fa

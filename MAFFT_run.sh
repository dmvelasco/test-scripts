#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fasta/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-mafft_run.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-mafft_run.txt
#SBATCH -J mafft
#SBATCH -a 2,5,9471,15800,21965,24852
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 20:00:00
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

# number of in array, with 1 CPU needed
# -a 1-27864%50
# why do I keep getting it wrong by using 27584?
#
# need up to 4 CPU
# -a 2,5,9471,15800,21965,24852

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare directories
dir1="/home/dmvelasc/bin"                                       # software binary directory
dir2="/group/jrigrp3/Velasco/Prunus/fasta/fasta-concat"         # multi-sequene FASTA  directory
dir3="/home/dmvelasc/Data/references/persica-SCF"               # FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/fasta/fasta-aligned"
dir5="/scratch/dmvelasc/fasta-aligned"

# basic set up
# input is multi-sequence FASTA

##### STEP 1: CREATE ARRARY OF GENE IDs #####
# creates array from gene ID file, each line (gene) a separate array item
# can use to create array job
# 27864 genes in gene ID list file

mkdir -p /scratch/dmvelasc/fasta-aligned

mapfile -s "$i" -n 1 -t gene < "$dir3"/Prunus_persica_v1.0_genes_list.gff3
#mapfile -t gene < "$dir3"/Prunus_persica_v1.0_genes_list.gff3
echo -e "${gene[0]}"

##### STEP 2: ALIGN MULTI FASTA SEQUENCE #####
#"$dir1"/mafft --localpair --maxiterate 1000 --inputorder "$dir2"/"${gene[0]}"_gene.fa > "$dir5"/"${gene[0]}"_gene_aln.fa
"$dir1"/mafft --localpair --maxiterate 1000 --inputorder "$dir2"/"${gene[0]}"_cds.fa > "$dir5"/"${gene[0]}"_cds_aln.fa

#mv "$dir5"/"${gene[0]}"_gene_aln.fa "$dir4"/
mv "$dir5"/"${gene[0]}"_cds_aln.fa "$dir4"/

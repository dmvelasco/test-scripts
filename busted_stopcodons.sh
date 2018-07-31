#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Script
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stopcodons-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stopcodons-stderr.txt
#SBATCH -J busted
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -a 1-20%2
#SBATCH --mem=8G
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of genes = number of arrays = 27864

### Load modules ###
module load hyphy

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
ref="/home/dmvelasc/Data/references/persica-SCF"			# reference directory
gene_pos_list="${ref}/Prunus_persica_v1.0_gene_position_list.txt"	# gene position list
FASTAdir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-aligned"		# final fasta directory
script="/home/dmvelasc/Projects/Prunus/Script"				# script directory that includes batch file and consensus tree
batch="/share/apps/hyphy-2.3.13/lib/TemplateBatchFiles"			# BUSTED batch file directory, trying to overcome misplacement of libv3 directory
finaldir="/home/dmvelasc/Projects/Prunus/Analysis/selection"
#### sample ID file
# column 1: ID, column2: other ID/information
#list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

echo -e "extract sample ID from Script/sample.txt using mapfile"
date

mapfile -s "$i" -n 1 -t line < "${gene_pos_list}"
# -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
# line is the array name (anything in this position is the array name)

# create array from line
locus=(`echo "${line[0]}"`)
# declare gene ID variable from array
gene_id="${locus[4]}"

# create SCRATCH DIRECTORY for temporary file placement
mkdir -p /scratch/dmvelasc/busted

#### Index BAM file
echo -e "Starting BUSTED analysis"
date

##### B U S T E D #####
#HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batch}"/CleanStopCodons.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln.fa" "${script}/splitstree_all.tree" "No/Yes" "" > "$finaldir"/"${gene_id}_cds_nostop.txt"
HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batch}"/CleanStopCodons.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln.fa" "No/Yes" "" > "$finaldir"/"${gene_id}_cds_nostop.txt"
# what is output directory? working directory? what is output file?

echo -e "BUSTED analysis finished"
date

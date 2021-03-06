#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Script
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stopcodons-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stopcodons-stderr.txt
#SBATCH -J busted
#SBATCH -p bigmemm
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -a 1-1000%5
#SBATCH --mem=8G
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of genes = number of arrays = 27864
# test with 1000

### Load modules ###
module load hyphy

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
ref="/home/dmvelasc/Data/references/persica-SCF"			# reference directory
gene_pos_list="${ref}/Prunus_persica_v1.0_gene_position_list.txt"	# gene position list
FASTAdir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-concat"		# final fasta directory
nostopdir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-nostop"
script="/home/dmvelasc/Projects/Prunus/Script"				# script directory that includes batch file and consensus tree
batch="/share/apps/hyphy-2.3.13/lib/TemplateBatchFiles"			# BUSTED batch file directory, trying to overcome misplacement of libv3 directory
scratch="/scratch/dmvelasc/fasta-stopcodon"
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
mkdir -p "$scratch"

#### Index BAM file
echo -e "Starting HyPhy stop codon cleaning for ${gene_id}"
date

##### B U S T E D #####
#HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batch}"/CleanStopCodons.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln.fa" "${script}/splitstree_all.tree" "No/Yes" "" > "$scratch"/"${gene_id}_cds_nostop.txt"
#HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batch}"/CleanStopCodons.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln.fa" "No/Yes" "" > "$scratch"/"${gene_id}_cds_nostop.txt"
#HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batch}"/CleanStopCodons.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln.fa" "No/No" "$scratch"/"${gene_id}_cds_aln_nostop.fa" > "$scratch"/"${gene_id}_cds_nostop.txt"
HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batch}"/CleanStopCodons.bf "Universal" "${FASTAdir}/${gene_id}_cds.fa" "No/No" "$scratch"/"${gene_id}_cds_nostop.fa" > "$scratch"/"${gene_id}_nostop.txt"
# Does not seem to have an output directory, but in the interactive command line version the last item is for the output file (/path/<filename>.fa)

mv "$scratch"/"${gene_id}_cds_nostop.fa" "${nostopdir}"
mv "$scratch"/"${gene_id}_cds_nostop.txt" "${nostopdir}"

echo -e "HyPhy stop codon removal finished"
date

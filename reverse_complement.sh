#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fasta/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-revcomp-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-revcomp-stderr.txt
#SBATCH -J revcomp
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -a 9,10,13-28,30-34,36-37,39,41-47,49-57,62-67%5
#SBATCH --mem=8G
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

### all arrays: 1-5,7-10,13-28,30-34,36-37,39,41-47,49-57,62-67%2

########## Extract FASTA for each GENE/CDS interval ##########
### Load modules ###
module load zlib

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
ref="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
scratch="/scratch/dmvelasc"				# Scratch directory
script="/home/dmvelasc/Projects/Prunus/Script"		# Script directory

#### sample ID file
# column 1: ID, column2: other ID/information
list="${script}/sample.txt"

# gene/cds intervals
gene_list="Prunus_persica_v1.0_genes_list.gff3"
gene_pos_list="Prunus_persica_v1.0_gene_position_list.txt"

##### mapfile to extract sample ID, each line is array item
mapfile -s "$i" -n 1 -t id < "${list}"
# -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
# id is the array name (anything in this position is the array name)

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
acc="${arr[0]}"
echo -e "$acc"

# create scratch directory for temporary file placement
mkdir -p "$scratch"/"$acc"

###############################################################################
# see possble solution at http://seqanswers.com/forums/showthread.php?t=50008 #
###############################################################################
echo -e "Reverse complement opposite strand FASTAs"
while read p; do
  # create an array from each line
  locus=(`echo "$p"`)
  # declare variables, created from array
  gene_id="${locus[4]}"
  strand="${locus[3]}"
  # reverse complement reverse strand gene sequence
  if [ "$strand" = '-' ]; then
    "$script"/DNA_reverse_complement.pl "$acc"/"$gene_id"_"$acc"_gene.fa > "$scratch"/"$acc"/"$gene_id"_"$acc"_temp.fa
    echo ">${acc}" > "$scratch"/"$acc"/"$gene_id"_"$acc"_gene.fa
    tail -n +2 "$scratch"/"$acc"/"$gene_id"_"$acc"_temp.fa | fold -w 60 - >> "$scratch"/"$acc"/"$gene_id"_"$acc"_gene.fa
    mv "$scratch"/"$acc"/"$gene_id"_"$acc"_gene.fa "$acc"/

    "$script"/DNA_reverse_complement.pl "$acc"/"$gene_id"_"$acc"_cds.fa > "$scratch"/"$acc"/"$gene_id"_"$acc"_temp.fa
    echo ">${acc}" > "$scratch"/"$acc"/"$gene_id"_"$acc"_cds.fa
    tail -n +2 "$scratch"/"$acc"/"$gene_id"_"$acc"_temp.fa | fold -w 60 - >> "$scratch"/"$acc"/"$gene_id"_"$acc"_cds.fa
    mv "$scratch"/"$acc"/"$gene_id"_"$acc"_cds.fa "$acc"/
    rm "$scratch"/"$acc"/"$gene_id"_"$acc"_temp.fa
  fi
done < "$ref"/"$gene_pos_list"

echo -e "end FASTA reverse complement process for CDS and genes"
date

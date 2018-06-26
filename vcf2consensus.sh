#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -a 1
#SBATCH --mem=16G
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

########## EXTRACT FASTA for each CDS/GENE ##########
# split on fasta header
#csplit -f "$prefix"fa_ -s file.fa '/>/' {*}

### Load modules ###
module load zlib

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
bin="/home/dmvelasc/bin"				# program directory
ref="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
scratch="/scratch/dmvelasc"				# Scratch directory
vcfcons="/home/dmvelasc/Software/vcftools/src/perl"	# VCFtools perl directory

# VCF file, sample IDS are PB01, PD02, etc.
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/all_jointcalls.vcf.gz"

# Asked for bgzip of VCF
#bgzip -c file.vcf > file.vcf.gz
#tabix -p vcf file.vcf.gz

#### sample ID file
# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

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

echo "begin conversion to FASTA"
date

##### Looks like it may be better to do this with samtools mpileup and then bcftools consensus

### Genes ###
echo -e "Create consensus FASTA for each gene region by:\n1. looping through list of gene coordinates and extracting full gene sequences, and then\n2. looping through CDS cds coordinates for each gene to create a single CDS FASTA."
date
###############################################################################
# see possble solution at http://seqanswers.com/forums/showthread.php?t=50008 #
###############################################################################
while read p; do
  # create an array from each line
  locus=(`echo "$p"`)
  # declare variables, created from array
  gene_interval="${locus[0]}:${locus[1]}-${locus[2]}"
  gene_id="${locus[4]}"
  # extract consensus gene sequence
  "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$gene_interval" | "$bin"/bcftools consensus "$vcf" -s "$acc" -o "$scratch"/"$acc"/"$gene_id"_"$acc".fa
done < "$ref"/"$gene_pos_list"

##### process each CDS from list of gene IDs and raw GATK alternate reference FASTA maker
echo -e "create consensus FASTA for CDS regions and concatenate to single FASTA file"
while read q; do
  touch "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa
  z=1
  # create CDS FASTA components from BAM
  while read r; do
    "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$r" | "$bin"/bcftools consensus "$vcf" -s "$acc" -o "$scratch"/"$acc"/"$q"_"$acc"_temp"$z".fa
    ((z++))
    tail -n +2 "$scratch"/"$acc"/"$q"_"$acc"_temp"$z".fa >> "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa
    rm "$scratch"/"$acc"/"$q"_"$acc"_temp"$z".fa
  done < "$ref"/cds_intervals/"$q".intervals

  # concatenate and create final CDS FASTA with basic file manipulations
  # concensus concatenation
  echo ">${acc}" > "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
  echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
  rm "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa
done < "$ref"/"$gene_list"

# move sample file directory from scratch
mv /scratch/dmvelasc/"$acc"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

echo "end CDS FASTA script"
date

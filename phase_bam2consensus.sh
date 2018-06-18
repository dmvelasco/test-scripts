#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -a 1-12%3
#SBATCH --mem=16G
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

########## PHASE BAMs & EXTRACT FASTA for each CDS/GENE ##########
# Part 1:
# 3. use samtools mpileup, bcftools call, bcftools consensus
# (initial: use samtools view and samtools fasta to extract each CDS/gene FASTA sequence,
# does not quite work as expected.)
# 4. append FASTA sequence to one file with appropriate FASTA header for each
# -------------------------------------------------
# Part 2: (separate)
# 1. create multi sequence alignment for each CDS/gene with mafft
# 2. selection analysis with BUSTED

### Load modules ###
module load zlib
module load bambam

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
bin="/home/dmvelasc/bin"				# program directory
ref="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
scratch="/scratch/dmvelasc"

# complete gff3: genes, CDS, UTR
genes="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_genes_only.gff3"
cds="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_cds_full.gff3"

#### sample ID file
# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

# gene/cds intervals
gene_list="Prunus_persica_v1.0_genes_list.gff3"
gene_pos_list="Prunus_persica_v1.0_gene_position_list.txt"

echo -e "extract sample ID from Script/sample.txt using mapfile"
date

mapfile -s "$i" -n 1 -t id < "${list}"
# -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
# id is the array name (anything in this position is the array name)

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
acc="${arr[0]}"
echo -e "$acc"

# create SCRATCH DIRECTORY for temporary file placement
mkdir -p "$scratch"/"$acc"_fasta/gene "$scratch"/"$acc"_fasta/cds

#### Index BAM file
echo -e "Check if HCrealign BAM file is indexed"
date

if [ ! -f "'$acc'_HCrealign.bam.bai" ]; then
  echo -e "Indexed file does not exist, indexing."
  "$bin"/samtools index "$acc"_HCrealign.bam
else
  echo -e "HCrealign BAM file index file exists, skipping."
fi

# Phase with SAMtools and index
echo "phase HC_realign BAM file with SAMtools"
#srun "$bin"/samtools calmd -AEur "$acc"_HCrealign.bam "$ref"/Prunus_persica_v1.0_scaffolds.fa | "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased
srun "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased "$acc"_HCrealign.bam
"$bin"/samtools index "$acc"_phased.0.bam
"$bin"/samtools index "$acc"_phased.1.bam


### Genes ###
echo -e "Process BAM file with bam2consensus to extract consensus FASTA file by GFF3 ID"
echo -e "Create consensus FASTA for each gene region by:\n1. splitting FASTA resulting from bam2consensus, and \n2. separating the whole gene sequences so as to not include in CDS file, and then \n3. identifying and concatenating CDS files for each gene to create a single CDS FASTA."
date

# Output phased consensus sequences from each phased BAM; direct to file otherwise defaults to stdout
while read g; do
  # Output full gene FASTA
  grep "$g" "$genes" > gene_"$g".gff3
  bam2consensus -g gene_"$g".gff3 "$acc"_phased.0.bam > "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.0_genes_gff3.fa
  bam2consensus -g gene_"$g".gff3 "$acc"_phased.1.bam > "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.1_genes_gff3.fa

  # Reformat full gene FASTA file
  # phase 0
  echo ">${g}_${acc}_phased.0_gene" > "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.0_gene.fa
  tail -n +2 "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.0_genes_gff3.fa  >> "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.0_gene.fa
  rm "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.0_genes_gff3.fa
  # phase 1
  echo ">${g}_${acc}_phased.1_gene" > "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.1_gene.fa
  tail -n +2 "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.1_genes_gff3.fa  >> "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.1_gene.fa
  rm "$scratch"/"$acc"_fasta/gene/"$g"_"$acc"_phased.1_genes_gff3.fa

  rm gene_"$g".gff3 # removes temporary GFF3 file with gene interval

  # Output CDS FASTA
  grep "$g" "$cds" > cds_"$g".gff3
  bam2consensus -g cds_"$g".gff3 "$acc"_phased.0.bam > "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.0_cds_gff3.fa
  bam2consensus -g cds_"$g".gff3 "$acc"_phased.1.bam > "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.1_cds_gff3.fa

  # Concatenate CDS FASTA files
  # split on fasta header (phase 0), -f is prefix, -s is file to split, '/>/' is regex on which to split, {*} is to split for all occurrences
  csplit -f "$g".0_cds_ -s "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.0_cds_gff3.fa '/>/' {*}
  numCDS=( `ls "$g".0_cds_* | wc -l` )
  echo ">${g}_${acc}_phased.0_cds" > "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.0_cds.fa
    for i in "$g".0_cds_*; do
      [ -f "$i" ] || break # loop breaks if no matching file
      tail -n +2 "$i" >> "$scratch"/"$acc"_temp_cds.fa
    done
  tr -d '\n' < "$scratch"/"$acc"_temp_cds.fa | fold -w 60 - > "$scratch"/"$acc"_folded_cds.fa
  cat "$scratch"/"$acc"_folded_cds.fa >> "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.0_cds.fa
  rm "$g".0_cds_* "$scratch"/"$acc"_temp_cds.fa "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.0_cds_gff3.fa

  # split on fasta header (phase 1)
  csplit -f "$g".1_cds_ -s "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.1_cds_gff3.fa '/>/' {*}
  numCDS=( `ls "$g".1_cds_* | wc -l` )
  echo ">${g}_${acc}_phased.1_cds" > "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.1_cds.fa
    for i in "$g".1_cds_*; do
      [ -f "$i" ] || break # loop breaks if no matching file
      tail -n +2 "$i" >> "$scratch"/"$acc"_temp_cds.fa
    done
  tr -d '\n' < "$scratch"/"$acc"_temp_cds.fa | fold -w 60 - > "$scratch"/"$acc"_folded_cds.fa
  cat "$scratch"/"$acc"_folded_cds.fa >> "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.1_cds.fa
  rm "$g".1_cds_* "$scratch"/"$acc"_temp_cds.fa "$scratch"/"$acc"_fasta/cds/"$g"_"$acc"_phased.1_cds_gff3.fa
  rm cds_"$g".gff3 # removes temporary GFF3 file with cds intervals
done < "$ref"/"$gene_list"

# move sample file directory from scratch, remove intermediate temporary files
mv "$scratch"/"$acc"_fasta/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/
rm "$scratch"/"$acc"_folded_cds.fa

echo "end bam2consensus script"
date

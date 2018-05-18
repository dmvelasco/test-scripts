#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

########## PHASE BAMs & EXTRACT FASTA for each CDS/GENE ##########
# In a loop?
# Part 1:
# 3. use samtools mpileup, bcftools call, bcftools consensus
# (initial: use samtools view and samtools fasta to extract each CDS/gene FASTA sequence,
# does not quite work as expected.)
# 4. append FASTA sequence to one file with appropriate FASTA header for each
# -------------------------------------------------
# Part 2: (separate)
# 1. create multi sequence alignment for each CDS/gene with mafft
# 2. selection analysis

### Load modules ###
module load zlib

### Declare directories ###
bin="/home/dmvelasc/bin"				# program directory
ref="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
scratch="/scratch/dmvelasc"

#### sample ID file
# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

# gene/cds intervals
gene_list="Prunus_persica_v1.0_genes_list.gff3"
gene_pos_list="Prunus_persica_v1.0_gene_position_list.txt"

echo "Begin phasing script"

# loop through each bam file
for i in {0..1}; do
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

  #### Index BAM file
  if [ ! -f "'$acc'_HCrealign.bam.bai" ]; then
    echo -e "Indexed file does not exist, indexing."
    "$bin"/samtools index "$acc"_HCrealign.bam
  else
    echo -e "Indexed file exists, skipping to phasing."
  fi
  ##### Phase with SAMtools #####
####### ->  srun "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased "$acc"_HCrealign.bam ##### keep to run full script
  # out is working/default directory

  #### Index phased BAM files
  if [ ! -f "'$acc'_phased.0.bam.bai" ]; then
    "$bin"/samtools index "$acc"_phased.0.bam
    "$bin"/samtools index "$acc"_phased.1.bam
  fi


  echo "begin conversion to FASTA"
  date

  ##### Extract FASTA from each BAM with SAMtools #####

##### Looks like it may be better to do this with samtools mpileup and then bcftools consensus

  ### Genes ###
  while read p; do
    # create an array from each line
    locus=(`echo "$p"`)
    # declare variables, created from array
    gene_interval="${locus[0]}:${locus[1]}-${locus[2]}"
    gene_id="${locus[4]}"

    # do for each phased bam of the sample
    "$bin"/samtools mpileup -u -r "$gene_interval" -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa -o "$scratch"/"$acc"/"$gene_id"_"$acc"_0.fa
    "$bin"/samtools mpileup -u -r "$gene_interval" -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa -o "$scratch"/"$acc"/"$gene_id"_"$acc"_1.fa
##    "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam "$gene_interval" | "$bin"/samtools fasta - > "$scratch"/"$acc"/"$gene_id"_"$acc"_0.fa
##    "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam "$gene_interval" | "$bin"/samtools fasta - > "$scratch"/"$acc"/"$gene_id"_"$acc"_1.fa
  done < "$ref"/"$gene_pos_list"

  ##### process each CDS from list of gene IDs and raw GATK alternate reference FASTA maker
  # will partly depend on how samtools outputs a fasta

  while read q; do
    # create phased CDS FASTA components from BAM
    # concatenate and create final CDS FASTA with basic file manipulations
    touch "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa
    touch "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa
##    while read r; do
      "$bin"/samtools mpileup -v -R "$ref"/cds_intervals/"$q".intervals -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa - -o "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa
      "$bin"/samtools mpileup -v -R "$ref"/cds_intervals/"$q".intervals -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa - -o "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa
##      "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam "$r" | "$bin"/samtools fasta - >> "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa
##      "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam "$r" | "$bin"/samtools fasta - >> "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa
##    done < "$ref"/cds_intervals/"$q".intervals
##    echo ">${acc}_0" > "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
##    echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
#    awk '{printf $0;}' "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc".fa
##    echo ">${acc}_1" >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
##    echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
#    awk '{printf $0;}' "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc".fa
#    rm "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
  done < "$ref"/"$gene_list"

  # move sample file directory from scratch
  mv /scratch/dmvelasc/"$acc"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

  echo "end CDS FASTA script"
  date
done

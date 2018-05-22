#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
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
  echo -e "Check if HCrealign BAM file is indexed"
  date

  if [ ! -f "'$acc'_HCrealign.bam.bai" ]; then
    echo -e "Indexed file does not exist, indexing."
    "$bin"/samtools index "$acc"_HCrealign.bam
  else
    echo -e "HCrealign BAM file index file exists, skipping to phasing."
  fi

  ##### Phase with SAMtools #####
  echo "Begin phasing"
  date

  srun "$bin"/samtools calmd -AEur "$acc"_HCrealign.bam "$ref"/Prunus_persica_v1.0_scaffolds.fa | "$bin"/samtools phase -A -A -b "$acc"_phased - ##### keep to run full script
  # previously srun "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased "$acc"_HCrealign.bam
  # out is working/default directory

  #### Index phased BAM files
  echo -e "Indexing phased BAM files"
  date
#  if [ ! -f "'$acc'_phased.0.bam.bai" ]; then
  "$bin"/samtools index "$acc"_phased.0.bam
  "$bin"/samtools index "$acc"_phased.1.bam
#  fi


  echo "begin conversion to FASTA"
  date

  ##### Extract FASTA from each BAM with SAMtools #####
  echo -e "pile up sequences and call VCF for first phased BAM file, phased_0"
  date
  # phase 0
  "$bin"/samtools mpileup -uABRx -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools call -c -O z -o "$acc"_phased.0.vcf
  "$bin"/bcftools index -f "$acc"_phased.0.vcf
  echo -e "pile up sequences and call VCF for second phased BAM file, phased_1"
  date
  # phase 1
  "$bin"/samtools mpileup -uABRx -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools call -c -O z -o "$acc"_phased.1.vcf
  "$bin"/bcftools index -f "$acc"_phased.1.vcf

##### Looks like it may be better to do this with samtools mpileup and then bcftools consensus

  ### Genes ###
  echo -e "loop through list of gene positions and extract full gene sequences and cds with samtools and bcftools consensus"
  date
  while read p; do
    # create an array from each line
    locus=(`echo "$p"`)
    # declare variables, created from array
    gene_interval="${locus[0]}:${locus[1]}-${locus[2]}"
    gene_id="${locus[4]}"
    # extract consensus gene sequence from phase 0
    "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$gene_interval" | "$bin"/bcftools consensus "$acc"_phased.0.vcf -o "$scratch"/"$acc"/"$gene_id"_"$acc"_0.fa
    # extract consensus gene sequence from phase 1
    "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$gene_interval" | "$bin"/bcftools consensus "$acc"_phased.1.vcf -o "$scratch"/"$acc"/"$gene_id"_"$acc"_1.fa

##    "$bin"/samtools mpileup -u -r "$gene_interval" -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa -o "$scratch"/"$acc"/"$gene_id"_"$acc"_0.fa
##    "$bin"/samtools mpileup -u -r "$gene_interval" -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa -o "$scratch"/"$acc"/"$gene_id"_"$acc"_1.fa
##    "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam "$gene_interval" | "$bin"/samtools fasta - > "$scratch"/"$acc"/"$gene_id"_"$acc"_0.fa
##    "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam "$gene_interval" | "$bin"/samtools fasta - > "$scratch"/"$acc"/"$gene_id"_"$acc"_1.fa
  done < "$ref"/"$gene_pos_list"

  ##### process each CDS from list of gene IDs and raw GATK alternate reference FASTA maker

  while read q; do
    # create phased CDS FASTA components from BAM
    # concatenate and create final CDS FASTA with basic file manipulations
    touch "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa
    touch "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa
    x=1
    while read r; do
      "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$r" | "$bin"/bcftools consensus "$acc"_phased.0.vcf -o "$scratch"/"$acc"/"$q"_"$acc"_temp"$x"_0.fa
      "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$r" | "$bin"/bcftools consensus "$acc"_phased.1.vcf -o "$scratch"/"$acc"/"$q"_"$acc"_temp"$x"_1.fa
##      "$bin"/samtools mpileup -v -R "$ref"/cds_intervals/"$q".intervals -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa - -o "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa
##      "$bin"/samtools mpileup -v -R "$ref"/cds_intervals/"$q".intervals -A -B -R -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam | "$bin"/bcftools consensus -f "$ref"/Prunus_persica_v1.0_scaffolds.fa - -o "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa
##      "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam "$r" | "$bin"/samtools fasta - >> "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa
##      "$bin"/samtools view -u -T "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam "$r" | "$bin"/samtools fasta - >> "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa
      ((x++))
      tail -n +2 "$scratch"/"$acc"/"$q"_"$acc"_temp"$x"_0.fa >> "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa
      tail -n +2 "$scratch"/"$acc"/"$q"_"$acc"_temp"$x"_1.fa >> "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa
      rm "$scratch"/"$acc"/"$q"_"$acc"_temp"$x"_0.fa "$scratch"/"$acc"/"$q"_"$acc"_temp"$x"_1.fa
    done < "$ref"/cds_intervals/"$q".intervals

    echo ">${acc}_0" >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
    echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
##  awk '{printf $0;}' "$scratch"/"$acc"/"$q"_"$acc"_cds_0.fa | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc".fa
    echo ">${acc}_1" >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
    echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
##  awk '{printf $0;}' "$scratch"/"$acc"/"$q"_"$acc"_cds_1.fa | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc".fa
    rm "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa
  done < "$ref"/"$gene_list"

  # move sample file directory from scratch
  mv /scratch/dmvelasc/"$acc"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

  echo "end CDS FASTA script"
  date
done

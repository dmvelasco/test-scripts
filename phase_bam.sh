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
scratch="/scratch/dmvelasc"

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

#### Index BAM file
echo -e "Check if HCrealign BAM file is indexed"
date

if [ ! -f "'$acc'_HCrealign.bam.bai" ]; then
  echo -e "Indexed file does not exist, indexing."
  "$bin"/samtools index "$acc"_HCrealign.bam
else
  echo -e "HCrealign BAM file index exists, skipping to phasing."
fi

##### Phase with SAMtools #####
echo "Begin phasing"
date

#### keep one of these to run full script, ONLY IF PHASING WITH SAMTOOLS FIRST
#srun "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased "$acc"_sorted_markdup.bam
#srun "$bin"/samtools calmd -AEur "$acc"_sorted_markdup.bam "$ref"/Prunus_persica_v1.0_scaffolds.fa | "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased
#srun "$bin"/samtools calmd -AEur "$acc"_HCrealign.bam "$ref"/Prunus_persica_v1.0_scaffolds.fa | "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased
# previously... retry
#srun "$bin"/samtools phase -A -A -Q 20 -b "$acc"_phased "$acc"_HCrealign.bam
# next option to run samtools view within loop to select region directly from BAM


#### Index phased BAM files
#echo -e "Indexing phased BAM files"
#date
#if [ ! -f "'$acc'_phased.0.bam.bai" ]; then
#"$bin"/samtools index "$acc"_phased.0.bam
#"$bin"/samtools index "$acc"_phased.1.bam
#fi


echo "begin conversion to FASTA"
date

##### Extract FASTA from each BAM with SAMtools #####
# phase 0
#echo -e "pile up sequences and call VCF for first phased BAM file, phased_0"
#date
#"$bin"/samtools mpileup -uARxE -q 30 -Q 20 -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.0.bam | "$bin"/bcftools call -m -O z -o "$acc"_phased.0.vcf
#"$bin"/bcftools index -f "$acc"_phased.0.vcf
# phase 1
#echo -e "pile up sequences and call VCF for second phased BAM file, phased_1"
#date
#"$bin"/samtools mpileup -uARxE -q 30 -Q 20 -f "$ref"/Prunus_persica_v1.0_scaffolds.fa "$acc"_phased.1.bam | "$bin"/bcftools call -m -O z -o "$acc"_phased.1.vcf
#"$bin"/bcftools index -f "$acc"_phased.1.vcf

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
  # extract consensus gene sequence from phase 0
  "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$gene_interval" | "$bin"/bcftools consensus "$vcf" -s "$acc" -o "$scratch"/"$acc"/"$gene_id"_"$acc".fa
#  # extract consensus gene sequence from phase 0
#  "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$gene_interval" | "$bin"/bcftools consensus "$acc"_phased.0.vcf -s "$acc"_phased.0.bam -o "$scratch"/"$acc"/"$gene_id"_"$acc"_0.fa
#  # extract consensus gene sequence from phase 1
#  "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$gene_interval" | "$bin"/bcftools consensus "$acc"_phased.1.vcf -s "$acc"_phased.1.bam -o "$scratch"/"$acc"/"$gene_id"_"$acc"_1.fa
done < "$ref"/"$gene_pos_list"

##### process each CDS from list of gene IDs and raw GATK alternate reference FASTA maker
echo -e "create consensus FASTA for CDS regions and concatenate to single FASTA file"
while read q; do
  touch "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa
#  touch "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa
#  touch "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa
  z=1
  # create phased CDS FASTA components from BAM
  while read r; do
    "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$r" | "$bin"/bcftools consensus "$vcf" -s "$acc" -o "$scratch"/"$acc"/"$q"_"$acc"_temp"$z".fa
#    "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$r" | "$bin"/bcftools consensus "$acc"_phased.0.vcf -s "$acc"_phased.0.bam -o "$scratch"/"$acc"/"$q"_"$acc"_temp"$z"_0.fa
#    "$bin"/samtools faidx "$ref"/Prunus_persica_v1.0_scaffolds.fa "$r" | "$bin"/bcftools consensus "$acc"_phased.1.vcf -s "$acc"_phased.1.bam -o "$scratch"/"$acc"/"$q"_"$acc"_temp"$z"_1.fa
    tail -n +2 "$scratch"/"$acc"/"$q"_"$acc"_temp"$z".fa >> "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa
#    tail -n +2 "$scratch"/"$acc"/"$q"_"$acc"_temp"$z"_0.fa >> "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa
#    tail -n +2 "$scratch"/"$acc"/"$q"_"$acc"_temp"$z"_1.fa >> "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa
    rm "$scratch"/"$acc"/"$q"_"$acc"_temp"$z".fa
#    rm "$scratch"/"$acc"/"$q"_"$acc"_temp"$z"_0.fa "$scratch"/"$acc"/"$q"_"$acc"_temp"$z"_1.fa
    ((z++))
  done < "$ref"/cds_intervals/"$q".intervals

  # concatenate and create final CDS FASTA with basic file manipulations
  # concensus concatenation
  echo ">${acc}" > "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
  echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
#  # phase 0
#  echo ">${acc}_0" >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
#  echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
#  # phase 1
#  echo ">${acc}_1" >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
#  echo $(cat "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa) | fold -w 60 - >> "$scratch"/"$acc"/"$q"_"$acc"_cds.fa
  rm "$scratch"/"$acc"/"$q"_"$acc"_cds_temp.fa
#  rm "$scratch"/"$acc"/"$q"_"$acc"_cds_0_temp.fa "$scratch"/"$acc"/"$q"_"$acc"_cds_1_temp.fa
done < "$ref"/"$gene_list"

# move sample file directory from scratch
mv /scratch/dmvelasc/"$acc"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

echo "end CDS FASTA script"
date

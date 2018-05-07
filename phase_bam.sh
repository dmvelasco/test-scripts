#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p serial
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=60G
set -e
set -u

########## PHASE BAMs & EXTRACT FASTA for each CDS/GENE ##########
# In a loop?
# Part 1:
# 1. prep prefix information to cycle through each BAM file
# (using the GATK produced HCrealign.bam files)
# 2. use samtools phase to output at least two new phased BAM files
# 3. use samtools view and samtools fasta to extract each CDS/gene FASTA sequence
# 4. append FASTA sequence to one file with appropriate FASTA header for each
# -------------------------------------------------
# Part 2:
# 1. create multi sequence alignment for each CDS/gene with mafft
# 2. selection analysis

### Load modules ###
module load zlib

### Declare directories ###
bin="/home/dmvelasc/bin"				# program directory
BAMdir="/group/jrigrp3/Velasco/Prunus/BAM"		# BAM directory
ref="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory
scratch="/scratch/dmvelasc"

#### sample ID file
# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

# gene intervals
intervals="/home/dmvelasc/Data/references/persica-SCF/Ppv1gene.intervals"

loop through each bam file
for i in {0..2}; do
  ##### mapfile to extract sample ID, each line is array item
  mapfile -s "$i" -n 1 -t id < "${list}"
  # -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
  # id is the array name (anything in this position is the array name)

  # create an array from each two column line
  arr=(`echo "${id[0]}"`)

  # declare variables, created from array
  acc="${arr[0]}"

####### Don't need?
#  # make directory to copy reference
#  mkdir -p "$acc"_temp

####### Don't need
#  # Index reference to temp directories
#  cp "$dir2"/Prunus_persica_v1.0_scaffolds.fa.* "$acc"_temp/
#  cp "$dir5"/Prunus_persica_v1.0_scaffolds.fa "$acc"_temp/

###### Don't need?
#  # create scratch directory for temporary file placement
#  mkdir -p /scratch/dmvelasc/"$acc"

  ##### Phase with SAMtools #####
  srun "$bin"/samtools phase -A -b "$acc"_phased -Q 20 "$BAMdir"/"$acc".HCrealign.bam
  # out is working/default directory

  echo "begin CDS FASTA script"
  date

  ##### Extract FASTA from each BAM with SAMtools #####

  ##### mapfile to extract cds/gene region
  mapfile -s "$i" -n 1 -t genes < "${intervals}"
  # -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
  # id is the array name (anything in this position is the array name)

  # create an array from each line
  cds=(`echo "${genes[0]}"`)
#  # declare variables, created from array
#  interval="${cds[0]}"

  # do for each phased bam of the sample
  samtools "$bin"/samtools view -u "$acc"_phased.0.bam "$cds" | "$bin"/samtools fasta -
  samtools "$bin"/samtools view -u "$acc"_phased.1.bam "$cds" | "$bin"/samtools fasta -
#  samtools "$bin"/samtools view -u "$acc"_phased.chimeric.bam "$cds" | "$bin"/samtools fasta
	concatenate fastas for each region

  # move phased BAMs to BAM directory?

  ##### process each CDS

  while read p; do
    # basic manipulations to create final CDS FASTA
    echo ">${acc}" > "$dir4"/"$acc"/"$p"_"$acc".fa
    awk '{printf $0;}' "$dir4"/"$acc"/"$p"_"$acc"_cds.fa | fold -w 60 - >> "$dir4"/"$acc"/"$p"_"$acc".fa
    rm "$scratch"/"$acc"/"$p"_"$acc"_cds.fa
  done < "$ref"/Prunus_persica_v1.0_genes_list.gff3

  # move sample file directory from scratch
  mv /scratch/dmvelasc/"$acc"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

  echo "end CDS FASTA script"
  date
done


##### IS BELOW NEEDED?? #####
#Simple UNIX code for splitting FASTA on '>'
#see https://www.biostars.org/p/2226/
#Solutions offered by users Biomonika (Noolean) and new

# csplit -f "$prefix"fa_ -s file.fa '/>/' {*}

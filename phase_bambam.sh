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
mkdir -p "$scratch"/"$acc"/gene "$scratch"/"$acc"/cds

#### Index BAM file
echo -e "Check if HCrealign BAM file is indexed"
date

if [ ! -f "'$acc'_HCrealign.bam.bai" ]; then
  echo -e "Indexed file does not exist, indexing."
  "$bin"/samtools index "$acc"_HCrealign.bam
else
  echo -e "HCrealign BAM file index file exists, skipping to phasing."
fi

### Genes ###
echo -e "Select region of BAM file and process with hapHunt"
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
  # phase and output FASTA
  srun "$bin"/samtools view -h -o "$gene_id"_"$acc"_gene.bam "$acc"_HCrealign.bam "$gene_interval"
  "$bin"/samtools index "$gene_id"_"$acc"_gene.bam
  hapHunt "$gene_id"_"$acc"_gene.bam > "$scratch"/"$acc"/gene/"$gene_id"_"$acc".fasta ### output to working directory, can this be redirected?
  rm "$gene_id"_"$acc"_gene.bam "$gene_id"_"$acc"_gene.bam.bai
done < "$ref"/"$gene_pos_list"

##### process each CDS from list of gene IDs
echo -e "create consensus FASTA for CDS regions and concatenate to single FASTA file"
while read q; do
  # create CDS FASTA components from BAM
  touch "$q"_"$acc"_bamlist.txt
  while read r; do
    srun "$bin"/samtools view -o "$q"_"$r"_"$acc"_cds.bam "$acc"_HCrealign.bam "$r"
    "$bin"/samtools index "$q"_"$r"_"$acc"_cds.bam
    echo -e ""$q"_"$r"_"$acc"_cds.bam" >> "$q"_"$acc"_bamlist.txt
  done < "$ref"/cds_intervals/"$q".intervals

  # concatenate BAM files for each CDS and perform hapHunt
  srun "$bin"/samtools cat -b "$q"_"$acc"_bamlist.txt -o "$q"_"$acc"_cds.bam
  "$bin"/samtools index "$q"_"$acc"_cds.bam
  hapHunt "$q"_"$acc"_cds.bam > "$scratch"/"$acc"/cds/"$gene_id"_"$acc".fasta

  # remove intermediate BAM files
  rm "$q"_*_"$acc"_cds.bam "$q"_*_"$acc"_cds.bam.bai

done < "$ref"/"$gene_list"

# move sample file directory from scratch
mv /scratch/dmvelasc/"$acc"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

echo "end CDS FASTA script"
date

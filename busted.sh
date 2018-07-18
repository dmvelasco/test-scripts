#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-phase-stderr.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -a 9
#SBATCH --mem=16G
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u




########## PHASE BAMs & EXTRACT FASTA for each CDS/GENE ##########
# Part 1:
# 3. use samtools mpileup, bcftools call, bcftools consensus
# (initial: use samtools view and samtools fasta to extract each CDS/gene FASTA sequence,
# does not quite work as expected.)
# 4. append FASTA sequence to one file with appropriate FASTA header for each sample

### Load modules ###
module load hyphy
module load zlib
module load bambam

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
bin="/home/dmvelasc/bin"					# program directory
ref="/home/dmvelasc/Data/references/persica-SCF"		# reference directory
final="/group/jrigrp3/Velasco/Prunus/fasta"			# final fasta directory
gene_pos_list="${ref}/Prunus_persica_v1.0_gene_position_list.txt"	# gene position list

#### sample ID file
# column 1: ID, column2: other ID/information
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

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
mkdir -p /home/dmvelasc/Projects/Prunus/Data/fasta/"$acc"_scaffolds

#### Index BAM file
echo -e "Check if ${acc}_HCrealign BAM file is indexed"
date



##########################################
############### B U S T E D ##############
##########################################

# EXAMPLE of putting BUSTED options on the command line with a batch file
$HYPHYMP /usr/local/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf "Universal" "/media/genomics/DATA1/AnophelesProteome/Transcriptome/OG1_5_3520_aln.fasta" "/media/genomics/DATA1/AnophelesProteome/Transcriptome/PAML/test.tree" "All" ""





if [ ! -f "${acc}_HCrealign.bam.bai" ]; then
  echo -e "Indexed file does not exist, indexing."
  "$bin"/samtools index "$acc"_HCrealign.bam
else
  echo -e "${acc}_HCrealign BAM file index file exists, skipping to phasing."
fi

mkdir -p "$final"/fasta-phased/"$acc"

echo -e "Extract phased FASTA for each gene ID (scaffold) with hapHunt"
date

##### do samtools view with header for each gene creating a new BAM file
##### then index and try bambam hapHunt to phase
# "$bin"/samtools view -bh -o temp_"$acc".bam "$acc"_HCrealign.bam "$gene_interval"

# total genese 27864
for i in {1..104}; do
  mapfile -s "$i" -n 1 -t line < "${gene_pos_list}"
  # -s number of rows to skip | -n number of rows to read | -t (remove leading/trailing whitespace?)
  # line is the array name (anything in this position is the array name)

  # create array from line
  locus=(`echo "${line[0]}"`)
  # declare variables, created from array
  gene_interval="${locus[0]}:${locus[1]}-${locus[2]}"
  gene_id="${locus[4]}"
  chr="${locus[0]}"

  # create a new BAM file with the selected scaffold
  srun "$bin"/samtools view -bh -o "$acc"_"$gene_id".bam "$acc"_HCrealign.bam "$gene_interval"
  # index the new BAM file
  "$bin"/samtools index "$acc"_"$gene_id".bam
  # phase the selected scaffold
  hapHunt "$acc"_"$gene_id".bam
  ### outputs to working directory, not easily redirected
  # select region of interest
  "$bin"/samtools faidx "$chr".fasta
  "$bin"/samtools faidx -c "$chr".fasta "${acc}_${gene_id}.bam.0:${locus[1]}-${locus[2]}" >> "${acc}_${gene_id}".fa # format like fasta_header:start-stop
  "$bin"/samtools faidx -c "$chr".fasta "${acc}_${gene_id}.bam.1:${locus[1]}-${locus[2]}" >> "${acc}_${gene_id}".fa # format like fasta_header:start-stop
# "If regions are specified, the subsequences will be retrieved and printed to stdout in the FASTA format"
# does not include FASTA headers for each region when multiple regions are specified
# does incude a header either if only one region specified or directed to stdout

  # move and rename file
  mv "$acc"_"$gene_id".fa "$final"/fasta-phased/"$acc"/
  # remove selected scaffold and associated index file
  rm "$acc"_"$gene_id".bam "$acc"_"$gene_id".bam.bai scaffold_*.fasta
done # < "$ref"/"$gene_pos_list"

echo "hapHunt FASTA extraction from BAM file finished"
date

# need to examine output and determine how to extract desired locus

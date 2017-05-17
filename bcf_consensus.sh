#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-vcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-vcf.txt
#SBATCH -J qry
#SBATCH -p serial
#SBATCH -a 1-7
#SBATCH -n 1
#SBATCH -c 1
set -e
set -u

x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Load zlib 1.2.8
module load zlib

# Declare BAM file list array
declare -a name=(PD03 PP04 PR01 PU01 PV02 PT01 PC01)

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"

# process sorted gff3 to acquire gene positions
# grep gene "$dir3"/Prunus_persica_v1.0_genes_sorted.gff3 | awk 'BEGIN{OFS="\t";} {print $1,$4,$5,$7,$9;}' - > "$dir3"/Prunus_persica_v1.0_gene_positions.txt
# process to split and get only gene ID in order to use as FASTA name
# awk 'BEGIN{OFS="\t";} {split($5,a,"[=.;]");} {print $1,$2,$3,$4,a[5];}' Prunus_persica_v1.0_gene_positions.txt > Prunus_persica_v1.0_gene_position_list.txt

##### BCFtools options for extracting information from vcf files #####
# consensus - each sample must be done separately
# create array of samples and loop across them?
"$dir1"/bcftools index -f "$dir2"/prunus-A.flt_test.vcf.bzip
"$dir1"/samtools faidx "$dir3"/Prunus_persica_v1.0_scaffolds.fa scaffold_5:11500000-11750000 | "$dir1"/bcftools consensus -i -s "${name["$i"]}".bam "$dir2"/prunus-A.flt_test.vcf.bzip -o "$dir2"/"${name["$i"]}"_test.fa

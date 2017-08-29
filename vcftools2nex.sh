#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-snapp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-snapp.txt
#SBATCH -J vcf2nex
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# Load zlib 1.2.8
module load zlib
module load java
module load vcftools

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"

# Declare other variables
infile="all_jointcalls.vcf"
prefix="test"
thin="5000"				# spacing between each SNP
matrix="split_alltest_final.txt"	# transposed file
outfile="split_alltest_final.nex"	# final file name

# NEXUS variables
#DIMENSIONS (see below)
#FORMAT
datatype="standard" #from Integer, which is what SNAPP uses
symbols='"012"'
interleave="no"
missing="?"

##### BCFtools options for extracting information from vcf files #####
# filter vcf for chromosomes 1-8, snps only

vcftools --vcf "$infile" --out "$prefix" \
--chr scaffold_1 \
--chr scaffold_2 \
--chr scaffold_3 \
--chr scaffold_4 \
--chr scaffold_5 \
--chr scaffold_6 \
--chr scaffold_7 \
--chr scaffold_8 \
--thin "$thin"\
--remove-indels \
--max-alleles 2 \
--min-meanDP 5 \
--minQ 250 \
--minGQ 30 \
--012

#This option (--012) outputs the genotypes as a large matrix.
#File ".012" contains the genotypes of each individual on a separate line (0, 1 and 2 represent the number of non-reference alleles, missing are -1).
#File ".012.indv" details the individuals in the main file.
#File ".012.pos" details the site locations in the main file.

########## "NEED TO REMOVE FIRST COLUMN FROM GENOTYPE (.012) FILE B/C JUST ROW #"
########## "REPLACE -1 WITH ?"
########## "CONCATENATE COLUMNS, CURRENTLY SEPARATED BY TAB"

# create modified .012 file
perl # replace -1 with ?, pipe to remove all tabs
perl -plne 's/-1/?/g' "$out1.vcf" | perl -plne 's/\t//g' - > "$out2".vcf

########## SPLIT NAMES IN ".012.indv" FILE AND ADD UNDERSCORE" ##############
########## CONCATENATE FILES BY LINE: paste file1.txt file2.txt

awk '{a=substr($1, 1, 2); b=substr($1, 3, 2); print a"_"b;}' "$prefix".012.indv | paste - "$prefix".012 > "$matrix"

# Nexus file variables requiring intermediate files
##### number of taxa, use .012.indv file or intermediate matrix file
Ntax=$(cat "$matrix" | wc -l)

##### number of snps, use modified .012 file to count characters then divide by # individuals
Nchar=$(( `cat "$prefix".012 | wc -c`/$Ntax )) #number of snps ------------ USE .012 FILE TO COUNT CHARACTERS THEN DIVIDE BY # INDIVIDUALS ----------

# Nexus file header and body creation
echo -e "#NEXUS\n[Written $(date)]\nBEGIN Data;" > "$outfile"
echo -e "\tDIMENSIONS NTAX="$Ntax" NCHAR="$Nchar";" >> "$outfile"
echo -e "\tFORMAT DATATYPE="$datatype" Symbols="$symbols" INTERLEAVE="$interleave" missing="$missing";" >> "$outfile"
echo "Matrix" >> "$outfile"
cat "$matrix" >> "$outfile"
echo -e ";\nEND;" >> "$outfile"

#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-hapcompass.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-hapcompass.txt
#SBATCH -J test
#SBATCH -p bigmemm
#SBATCH -t 14-00:00:00
#SBATCH -a 1-2
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=32G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft and WhatsHap)
module load java
module load vcftools
module load tabix

#############################
### Set up the parameters ###
#############################

####### PATHS #######
# HapCompass program
hapcompass="/home/dmvelasc/Software/hapcompass/hapcompass.jar"
# genome reference file
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# Joint VCF file - test joint VCF file
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/all_jointcalls.vcf"
# filtered joint VCF file - dulcis test VCF file
vcf_dir="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# SMC++ prepped file directory
bam_dir="/group/jrigrp3/Velasco/Prunus/BAM"

####################
### Begin script ###
####################

##### CREATE SAMPLE PREFIX #####
# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# path to sample list
list="/home/dmvelasc/Projects/Prunus/Script/sample_vcf.txt"

# mapfile to extract sample ID and read name information, each line is array item
mapfile -s "$i" -n 1 -t id < "${list}"
# -s number of rows to skip
# -n number of rows to read
# -t (cannot remember, but does not correspond to last item; remove leading/trailing whitespace?)
# id is the array name (anything in this position is the array name if nothing then called ?array)

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
sample="${arr[0]}"

# Can only do one genotype at a time, need to filter with VCFtools for each genotype
vcftools --vcf "$vcf" --out "$sample" --indv "$sample" --recode

# HapCompass script, export hard java "heap memory" limit range
export _JAVA_OPTIONS="-Xms4g -Xmx30g"
#java -Xmx30g -jar "$hapcompass" -bam "$bam_dir"/"$sample"_sorted_markdup.bam -vcf "$vcf_dir"/"$sample".recode.vcf -o "$vcf_dir"/"$sample"_phased.vcf
java -Xmx30g -jar "$hapcompass" -bam "$bam_dir"/"$sample"_HCrealign.bam -vcf "$vcf_dir"/"$sample".recode.vcf -o "$vcf_dir"/"$sample"_phased.vcf

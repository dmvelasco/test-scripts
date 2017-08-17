#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-GATK-eval.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-GATK-eval.txt
#SBATCH -p bigmemh
#SBATCH -a 55
#SBATCH -J eval
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=24000

set -e
set -u

# running to produce an index file for PD10
# 54 should be PD10

module load samtools/1.3.1
module load bamtools
module load java/1.8
# load GATK dependencies
module load R/3.3.1
module load maven/3.2.3
#module load GATK/3.6

########## WRITTEN BY D. VELASCO ###########

##########################################################################################################################
### picard verion: 2.9                                                                                                 ###
### GATK version: 3.7                                                                                                  ###
##########################################################################################################################

#############################
### Set up the parameters ###
#############################
# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# number of threads
#threads=6
# location for the picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location for the GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"
# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# genome reference index file location
gindex="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa.fai"
# genome dictionary file location
dictionary="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.dict"
# VCF directory
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# BAM directory
BAM="/group/jrigrp3/Velasco/Prunus/BAM"
# scratch directory for temporary storage
scratch="/scratch/dmvelasc"

##### CREATE SAMPLE PREFIX #####
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


###########
######### mk sample subdirectory and copy reference?

echo "########## evaluate GVCF ##########";
date

java -Xmx20g -jar "$GATK" -T VariantEval \
    -R "$genome" \
    -o "$sample".eval.grp \
    --eval:set1 "$sample".g.vcf

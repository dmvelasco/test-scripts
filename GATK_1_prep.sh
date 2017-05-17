#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-GATK1prep.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-GATK1prep.txt
#SBATCH -p serial
#SBATCH -J GATK1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 4-00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=3000M

set -e
set -u

module load samtools/1.3.1
module load bamtools
module load java/1.8
# load GATK dependencies
module load R/3.3.1
module load maven/3.2.3
#module load GATK/3.6

#######################################################################################
### picard verion: 2.9.0                                                            ###
### GATK version: 3.7                                                               ###
#######################################################################################

echo "########## set up directories  ##########";
date

###############################
###   Set up directories    ###
###############################
bwa="/home/dmvelasc/bin/bwa"

echo "########## set up parameters  ##########";
date

###############################
###  Set up the parameters  ###
###############################

# location of picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location of GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"
# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"

###############################
### Step 1: Index reference ###
###############################

# Step1: BWA index for reference genome

echo "########## prepare the BWA index and picard/GATK dictionary files ##########";
date

echo "indexing reference genome with BWA";
date

#"$bwa" index -a bwtsw "$genome" #check that options are correct, however, already done so skip step

echo "creating picard/GATK dictionary";
date

java -Xmx3g -jar "$picard" CreateSequenceDictionary \
 R="$genome" \
 O=/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.dict

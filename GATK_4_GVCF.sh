#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-GATK-gVCF.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-GATK-gVCF.txt
#SBATCH -p bigmemm
#SBATCH -a 1-20%
#SBATCH -J GATK
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 10-00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=32G

set -e
set -u

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

# number of threads
#threads=6
# location for the picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location for the GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"
# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))
# full set of initial
declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# set that had different quality encoding
#declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01)
sample="${id["$i"]}"

# Declare prefix array
# Public sequences for later
#declare -a pub=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES


# Step1: validate file

echo "########## Validate BAM File ##########";
date

java -Xmx20g -jar "$picard" ValidateSamFile \
    I="$sample".bam \
    MODE=SUMMARY


## Step2: call variants
echo "########## call variants and output GVCF ##########";
date

# run for individual bams, combine in joint genotyping later
java -Xmx20g -jar "$GATK" -T HaplotypeCaller \
    -R "$genome" \
    -I "$sample".bam \
    -o "$sample".g.vcf \
    -bamout "$sample"_HCrealign.bam \
    -ERC GVCF
# Every time there is an update on the project simply re-run the quick GenotypeGVCFs step
# on all the samples available. The expensive and time-consuming part of calculating
# genotype likelihoods is thus only done once on each sample which reduces wasted time
# and resources for each additional sample added to the project.
    # -bamout option "allows you to ask for the realigned version of the bam"

#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-GATK-gVCF.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-GATK-gVCF.txt
#SBATCH -p serial
#SBATCH -a 1-9
#SBATCH -J GATK
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 5-00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=24000

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
#declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210)
#declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
#declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502985 SRR502994 SRR502992 SRR502990 SRR502987 SRR502986 SRR503000 SRR502983 SRR502997 SRR502983 SRR502997 SRR502995 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
#declare -a pub=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES


# Step1: validate file

#echo "########## Validate BAM File ##########";
#date

java -Xmx20g -jar "$picard" ValidateSamFile \
    I="$sample"_sorted_markdup.bam \
    MODE=SUMMARY


# HaplotypeCaller does local reassembly, below (with -bamout option) allows visualization or reassembly
# Step2: realign messy regions

#echo "########## Realign Messy Regions of the Genome ##########";
#date

#java -Xmx20g -jar "$GATK" -T HaplotypeCaller \
#    -R "$genome" \
#    -I "$sample"_realigned.bam \
#    -o "$sample"_debug.vcf \
#    -bamout "$sample"_out.bam \
#    -forceActive -disableOptimizations
    # -bamout option "allows you to ask for the realigned version of the bam"


# Step3: call variants

echo "########## call variants and output GVCF ##########";
date
# run for individual bams, combine in joint genotyping later

java -Xmx20g -jar "$GATK" -T HaplotypeCaller \
    -R "$genome" \
    -I "$sample"_sorted_markdup.bam \
    -o "$sample".g.vcf \
    -ERC GVCF
# Every time there is an update on the project simply re-run the quick GenotypeGVCFs step
# on all the samples available. The expensive and time-consuming part of calculating
# genotype likelihoods is thus only done once on each sample which reduces wasted time
# and resources for each additional sample added to the project.

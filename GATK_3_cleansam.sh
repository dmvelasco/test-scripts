#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-GATK3.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-GATK3.txt
#SBATCH -p med
#SBATCH -a 12,38
#SBATCH -J GATK
#SBATCH -n 1
#SBATCH -c 10
#SBATCH -t 7-00:00:00
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

########## D. VELASCO ###########

##################################################################################################################################
### picard verion: 2.9.0                                                                                                       ###
### GATK version: 3.7                                                                                                          ###
### bam --> add group name --> fixmate --> remove duplicate --> Realignment Around Indels --> Base Quality Score Recalibration ###
##################################################################################################################################

echo "########## begin the pipeline ##########";
date

#############################
### Set up the parameters ###
#############################
# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# location for the picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location for the GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"
# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# BAM directory
BAM="/group/jrigrp3/Velasco/Prunus/BAM"
# scratch directory for temporary storage
scratch="/scratch/dmvelasc"

##### CREATE SAMPLE PREFIX #####
# path to sample list
list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

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

# make scratch directory for job
mkdir -p /scratch/dmvelasc

####################################################################################################################################
### previous BWA-mem align(bwa mem -M) and get the bam file, and now we need to filter the bam file and get the clean alignment  ###
### The -M flag causes BWA to mark shorter split hits as secondary (essential for Picard compatibility)                          ###
####################################################################################################################################

# prior to running script
#samtools view PG01.bam | grep -v "HWI-EAS373:6:70:857:521" > PG01_fix.sam
#samtools view PP10.bam | grep -v "HWI-EAS373:6:63:513:758" > PP10_fix.sam

# Step1: add group name, GATK will not work without a read group tag ( it looks like that : @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 )
echo "Step 1: add RG group name and sort";
date

java -Xmx20g -jar $picard AddOrReplaceReadGroups \
 I="$sample"_fix.sam \
 O="$scratch"/"$sample"_with_RG.bam \
 SORT_ORDER=coordinate \
 RGID=$sample \
 RGLB=$sample \
 RGPL=illumina \
 RGPU=$sample \
 RGSM=$sample \
 CREATE_INDEX=true

mv "$scratch"/"$sample"_with_RG.* "$BAM"/

# Step2: remove duplicate reads
echo "Step 2: remove duplicate reads";
date

java -Xmx20g -jar $picard MarkDuplicates \
 VALIDATION_STRINGENCY=LENIENT \
 I="$sample"_with_RG.bam \
 O="$scratch"/"$sample"_sorted_markdup.bam \
 METRICS_FILE="$sample"_metrics.txt \
 CREATE_INDEX=true

mv "$scratch"/"$sample"_sorted_markdup.* "$BAM"/
rm "$sample"_fix.sam
rm "$sample"_with_RG.bam

# Step3: Realignment Around Indels
# indel realignment not required with GATK version 3.6 and later
# HaplotypeCaller performs indel realignment, see https://software.broadinstitute.org/gatk/blog?id=7847


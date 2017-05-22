#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-GATK3.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-GATK3.txt
#SBATCH -p bigmemm
#SBATCH -a 1-10%2
#SBATCH -J GATK
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 4-00:00
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

# location for the picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location for the GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))
#declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
declare -a id=(PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10 PS02)
sample="${id["$i"]}"

# Public sequences for later
#declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502985 SRR502994 SRR502992 SRR502990 SRR502987 SRR502986 SRR503000 SRR502983 SRR502997 SRR502983 SRR502997 SRR502995 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
#declare -a pub=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES

# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"

####################################################################################################################################
### previous BWA-mem align(bwa mem -M) and get the bam file, and now we need to filter the bam file and get the clean alignment  ###
### The -M flag causes BWA to mark shorter split hits as secondary (essential for Picard compatibility)                          ###
####################################################################################################################################

# Step1: add group name, GATK will not work without a read group tag ( it looks like that : @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 )
echo "Step 1: add RG group name and sort";
date

java -Xmx20g -jar $picard AddOrReplaceReadGroups \
 I="$sample".bam \
 O="$sample"_with_RG.bam \
 SORT_ORDER=coordinate \
 RGID=$sample \
 RGLB=$sample \
 RGPL=illumina \
 RGPU=$sample \
 RGSM=$sample \
 CREATE_INDEX=true

# Step2: remove duplicate reads
echo "Step 2: remove duplicate reads";
date

java -Xmx20g -jar $picard MarkDuplicates \
 VALIDATION_STRINGENCY=LENIENT \
 I="$sample"_with_RG.bam \
 O="$sample"_sorted_markdup.bam \
 METRICS_FILE="$sample"_metrics.txt \
 CREATE_INDEX=true

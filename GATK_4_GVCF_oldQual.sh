#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-GATK-gVCF.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-GATK-gVCF.txt
#SBATCH -p bigmemm
#SBATCH -a 56-67%10
#SBATCH -J GATK
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 10-00:00:00
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
mkdir -p "$scratch"/"$sample"
cp "$genome" "$scratch"/"$sample"/
cp "$dictionary" "$scratch"/"$sample"/

# Step1: validate file
echo "########## Index BAM File ##########";
date

# Check if BAM is indexed, if not then index and move to same location as BAM
file="${sample}_sorted_markdup.bai"
if [ -e "$file" ]
then
#    echo "$file found, skipping index build"
    :
else
#    echo "$file not found, building index"
    java -jar "$picard" BuildBamIndex \
        I="$sample"_sorted_markdup.bam
#    mv "$sample"_sorted_markdup.bam.bai final/ #apparently creates in same directory as BAM file
fi

echo "########## Validate BAM File ##########";
date

java -Xmx20g -jar "$picard" ValidateSamFile \
    I="$sample"_sorted_markdup.bam \
    MODE=SUMMARY

# Call variants
## Step2: call variants
echo "########## call variants and output GVCF ##########";
date

# run for individual bams, combine in joint genotyping later
java -Xmx20g -jar "$GATK" -T HaplotypeCaller \
    -R "$scratch"/"$sample"/Prunus_persica_v1.0_scaffolds.fa \
    -fixMisencodedQuals \
    -I "$sample"_sorted_markdup.bam \
    -o "$scratch"/"$sample".g.vcf \
    -bamout "$scratch"/"$sample"_HCrealign.bam \
    -ERC GVCF
# Every time there is an update on the project simply re-run the quick GenotypeGVCFs step
# on all the samples available. The expensive and time-consuming part of calculating
# genotype likelihoods is thus only done once on each sample which reduces wasted time
# and resources for each additional sample added to the project.
    # -bamout option "allows you to ask for the realigned version of the bam"

##### Step 3: cleanup
rm "$scratch"/"$sample"/Prunus_persica_v1.0_scaffolds.fa
rm "$scratch"/"$sample"/Prunus_persica_v1.0_scaffolds.dict"
mv "$scratch"/"$sample"_HCrealign.bam /group/jrigrp3/Velasco/Prunus/BAM/
mv "$scratch"/"$sample".* "$vcf"/
rmdir "$scratch"/"$sample"

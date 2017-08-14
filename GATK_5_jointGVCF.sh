#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-GATK_jointgvcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-GATK_jointgvcf.txt
#SBATCH -p bigmemm
#SBATCH -J GATK5
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 4-00:00
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

########################################################################################################
### picard verion: 2.9                                                                               ###
### GATK version: 3.7                                                                                ###
########################################################################################################

#############################
### Set up the parameters ###
#############################

#threads
#threads=16
# location for the picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location for the GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"
# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"


# Declare sample ID numbers/prefix array
declare -a id=(PB01 PR01 PU01 PC01 PV01 PV02 PV03 PV04 PV05 PV06 PF01 PG02 PG03 PG04 PG05 PS01 PS02 PS03 PS04 PK01 PM01 PM02 PM03 PM04 PM05 PM06 PT01 PP01 PP02 PP03 PP04 PP05 PP06 PP08 PP09 PP11 PP12 PP13 PP14 PP15 PP37 PP38 PP39 PP40 PD01 PD02 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10 PD11 PD12 PD13 PD14 PD16 PD17 PD18 PD19 PD20 PD21)
# missing PG01, PP07, PP10
# PD10 not indexed

#declare -a id=(PC01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PR01 PS02) # test set

####################################################################################################################################

# joint genotyping using GVCF files

echo "########## Joint Genotyping with GVCFs ##########";
date

# test run
#java -Xmx20g -jar "$GATK" -T GenotypeGVCFs \
#    -R "$genome" \
#    -V "${id[0]}".g.vcf \
#    -V "${id[1]}".g.vcf \
#    -V "${id[2]}".g.vcf \
#    -V "${id[3]}".g.vcf \
#    -V "${id[4]}".g.vcf \
#    -V "${id[5]}".g.vcf \
#    -V "${id[6]}".g.vcf \
#    -V "${id[7]}".g.vcf \
#    -V "${id[8]}".g.vcf \
#    -V "${id[9]}".g.vcf \
#    -o test_jointcalls.vcf

java -Xmx20g -jar "$GATK" -T GenotypeGVCFs \
    -R "$genome" \
    -V "${id[0]}".g.vcf \
    -V "${id[1]}".g.vcf \
    -V "${id[2]}".g.vcf \
    -V "${id[3]}".g.vcf \
    -V "${id[4]}".g.vcf \
    -V "${id[5]}".g.vcf \
    -V "${id[6]}".g.vcf \
    -V "${id[7]}".g.vcf \
    -V "${id[8]}".g.vcf \
    -V "${id[9]}".g.vcf \
    -V "${id[10]}".g.vcf \
    -V "${id[11]}".g.vcf \
    -V "${id[12]}".g.vcf \
    -V "${id[13]}".g.vcf \
    -V "${id[14]}".g.vcf \
    -V "${id[15]}".g.vcf \
    -V "${id[16]}".g.vcf \
    -V "${id[17]}".g.vcf \
    -V "${id[18]}".g.vcf \
    -V "${id[19]}".g.vcf \
    -V "${id[20]}".g.vcf \
    -V "${id[21]}".g.vcf \
    -V "${id[22]}".g.vcf \
    -V "${id[23]}".g.vcf \
    -V "${id[24]}".g.vcf \
    -V "${id[25]}".g.vcf \
    -V "${id[26]}".g.vcf \
    -V "${id[27]}".g.vcf \
    -V "${id[28]}".g.vcf \
    -V "${id[29]}".g.vcf \
    -V "${id[30]}".g.vcf \
    -V "${id[31]}".g.vcf \
    -V "${id[32]}".g.vcf \
    -V "${id[33]}".g.vcf \
    -V "${id[34]}".g.vcf \
    -V "${id[35]}".g.vcf \
    -V "${id[36]}".g.vcf \
    -V "${id[37]}".g.vcf \
    -V "${id[38]}".g.vcf \
    -V "${id[39]}".g.vcf \
    -V "${id[40]}".g.vcf \
    -V "${id[41]}".g.vcf \
    -V "${id[42]}".g.vcf \
    -V "${id[43]}".g.vcf \
    -V "${id[44]}".g.vcf \
    -V "${id[45]}".g.vcf \
    -V "${id[46]}".g.vcf \
    -V "${id[47]}".g.vcf \
    -V "${id[48]}".g.vcf \
    -V "${id[49]}".g.vcf \
    -V "${id[50]}".g.vcf \
    -V "${id[51]}".g.vcf \
    -V "${id[52]}".g.vcf \
    -V "${id[53]}".g.vcf \
    -V "${id[54]}".g.vcf \
    -V "${id[55]}".g.vcf \
    -V "${id[56]}".g.vcf \
    -V "${id[57]}".g.vcf \
    -V "${id[58]}".g.vcf \
    -V "${id[59]}".g.vcf \
    -V "${id[60]}".g.vcf \
    -V "${id[62]}".g.vcf \
    -V "${id[63]}".g.vcf \
   -o all_jointcalls.vcf

# Every time there is an update on the project simply re-run the quick GenotypeGVCFs step
# on all the samples available. The expensive and time-consuming part of calculating
# genotype likelihoods is thus only done once on each sample which reduces wasted time
# and resources for each additional sample added to the project.

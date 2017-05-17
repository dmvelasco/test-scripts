#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-GATK_jointgvcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-GATK_jointgvcf.txt
#SBATCH -p bigmemh
#SBATCH -J GATK
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
#declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
declare -a id=(PS02 PD01 PP15)

# Public sequences for later - REMEMBER SOME HAVE MULTIPLE SRA RUNS
#declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502999 SRR502985 SRR502994 SRR502992 SRR502993 SRR502990 SRR502991 SRR502987 SRR502989 SRR502986 SRR503000 SRR503001 SRR502983 SRR502997 SRR502995 SRR502996 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
#declare -a pub=(PD11 PD12 PD13 PD14 PG01.1 PG01.2 PP01 PP02 PP03.1 PP03.2 PP04.1 PP04.2 PP05.1 PP05.2 PP06 PP07.1 PP07.2 PP08 PP09 PP10.1 PP10.2 PP11 PP12 PP13 PP14 PS01 PV01)
#declare -a id=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES

####################################################################################################################################

# joint genotyping using GVCF files

echo "########## Joint Genotyping with GVCFs ##########";
date

java -Xmx20g -jar "$GATK" -T GenotypeGVCFs \
    -R "$genome" \
    -V "${id[0]}".g.vcf \
    -V "${id[1]}".g.vcf \
    -V "${id[2]}".g.vcf \
    -o test_jointcalls.vcf

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
#    -V "${id[10]}".g.vcf \
#    -V "${id[11]}".g.vcf \
#    -V "${id[12]}".g.vcf \
#    -V "${id[13]}".g.vcf \
#    -V "${id[14]}".g.vcf \
#    -V "${id[15]}".g.vcf \
#    -V "${id[16]}".g.vcf \
#    -V "${id[17]}".g.vcf \
#    -V "${id[18]}".g.vcf \
#    -V "${id[19]}".g.vcf \
#    -o all_jointcalls.vcf

# Every time there is an update on the project simply re-run the quick GenotypeGVCFs step
# on all the samples available. The expensive and time-consuming part of calculating
# genotype likelihoods is thus only done once on each sample which reduces wasted time
# and resources for each additional sample added to the project.

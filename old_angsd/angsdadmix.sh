#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-2pop-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-2pop-stderr.txt
#SBATCH -J admix
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# Array variables
#x=$SLURM_ARRAY_TASK_ID

# Load zlib 1.2.8
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"				# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Script"		# script directory
dir3="/home/dmvelasc/Data/references/persica-SCF"	# reference genome FASTA directory
dir4="/home/dmvelasc/Data/references/cerasifera"	# ancestral genome FASTA directory
dir5="/home/dmvelasc/Software/angsd/misc"		# ANGSD subprogram directory
dir6="/home/dmvelasc/Software/ngsTools/ngsPopGen"	# ngsPopGen directory, subdirectory of ngsTools

#Declare other variables
taxa1="PD"				# file prefix for dulcis
taxa2="PP"				# file prefix for persica
taxa3="PAadmixtest2"			# file prefix for combined dulcis and persica
bams1="angsd.dulcisbam.txt"		# list of BAM files for dulcis
bams2="angsd.persicabam.txt"		# list of BAM files for persica
bams3="angsd.admixtestbam.txt"		# list of BAM files for dulcis and persica combined
ref="Prunus_persica_v1.0_scaffolds.fa"	# file name for reference genome
anc="PC01.fa"				# file name for ancestral genome
threads="8"				# number of threads used
minMapQ="30"				# minimum mapping quality
minQ="20"				# mimimum base quality
r="scaffold"				# region prefix

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

##### VARIABLES - don't think needed for admixture

#date
#echo "Prepping variables..."
#N_IND1=$(cat "$dir2"/"$bams1" | wc -l)
#N_IND2=$(cat "$dir2"/"$bams2" | wc -l)
#N_IND3=$(cat "$dir2"/"$bams3" | wc -l)
#echo -e "Number of $taxa1 individuals: $N_IND1 \nNumber of $taxa2 individuals: $N_IND2 \nNumber of $taxa3 individuals: $N_IND3"
#CHR1=$(( $N_IND1 * 2 ))
#CHR2=$(( $N_IND2 * 2 ))
#CHR3=$(( $N_IND3 * 2 ))
#echo -e "Number of $taxa1 chromosomes: $CHR1 \nNumber of $taxa2 chromosomes: $CHR2 \nNumber of $taxa3 chromosomes: $CHR3"
#date

#######################################
# COMBINED DULCIS AND PERSICA SAMPLES #
#######################################

##### ADMIXTURE - WORKS
date
echo "Calculating admixture..."
#example for creating Beagle input file, -doGlf 2; use combined bam file list
"$dir1"/angsd -GL 1 -out "$taxa3".admix -nThreads "$threads" -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam "$dir2"/"$bams3"
"$dir1"/NGSadmix -likes "$taxa3".admix.beagle.gz -K 2 -P "$threads" -o "$taxa3"_admix -minMaf 0.05
echo "Done calculating admixture."
date
########

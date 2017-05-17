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
taxa="PA"				# file prefix for combined dulcis and persica
bams="angsd.admix.txt"			# list of BAM files for dulcis and persica combined
ref="Prunus_persica_v1.0_scaffolds.fa"	# file name for reference genome
anc="PC01.fa"				# file name for ancestral genome
threads="8"				# number of threads used
minMapQ="30"				# minimum mapping quality
minQ="20"				# mimimum base quality
r="scaffold"				# region prefix

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

#######################################
# COMBINED DULCIS AND PERSICA SAMPLES #
#######################################

##### ADMIXTURE - WORKS
date
echo "Calculating admixture..."
#example for creating Beagle input file, -doGlf 2; use combined bam file list
#"$dir1"/angsd -GL 1 -out "$taxa".admix -nThreads "$threads" -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam "$dir2"/"$bams"
#"$dir1"/NGSadmix -likes "$taxa".admix.beagle.gz -K 2 -P "$threads" -o "$taxa"_admix -minMaf 0.05
"$dir1"/NGSadmix -likes "$taxa".admix.beagle.gz -K 3 -P "$threads" -o "$taxa"_admix_K3 -minMaf 0.05
"$dir1"/NGSadmix -likes "$taxa".admix.beagle.gz -K 4 -P "$threads" -o "$taxa"_admix_K4 -minMaf 0.05
"$dir1"/NGSadmix -likes "$taxa".admix.beagle.gz -K 5 -P "$threads" -o "$taxa"_admix_K5 -minMaf 0.05
"$dir1"/NGSadmix -likes "$taxa".admix.beagle.gz -K 6 -P "$threads" -o "$taxa"_admix_K6 -minMaf 0.05
echo "Done calculating admixture."
date
########

# WITH INBREEDING
#date
#echo "Calculating admixture with inbreeding..."
#example for creating Beagle input file, -doGlf 2; use combined bam file list
#"$dir1"/angsd -GL 1 -out "$taxa"_F.admix -nThreads "$threads" -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -indF "$taxa"_F.indF -bam "$dir2"/"$bams"
#"$dir1"/NGSadmix -likes "$taxa"_F.admix.beagle.gz -K 2 -P "$threads" -o "$taxa"_F_admix -minMaf 0.05
#echo "Done calculating admixture with inbreeding."
#date

#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A-prep2pop-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A-prep2pop-stderr.txt
#SBATCH -J prep2pop
#SBATCH -p bigmemm
#SBATCH -a 1
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Array variables
x=$SLURM_ARRAY_TASK_ID

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
taxa3="PE"				# file prefix for combined dulcis and persica
bams1="angsd.dulcisbam.txt"		# list of BAM files for dulcis
bams2="angsd.persicabam.txt"		# list of BAM files for persica
bams3="angsd.domesticatesbam.txt"	# list of BAM files for dulcis and persica combined
ref="Prunus_persica_v1.0_scaffolds.fa"	# file name for reference genome
anc="PC01.fa"				# file name for ancestral genome
threads="4"				# number of threads used
minMapQ="30"				# minimum mapping quality
minQ="20"				# mimimum base quality
r="scaffold"_"$x":			# region prefix

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

##### CALCULATING GENOTYPES, GENOTYPE LIKELIHOODS, SFS, MAFS, ETC. WHILE ACCOUNTING FOR INBREEDING

date
echo "Prepping variables..."
N_IND1=$(cat "$dir2"/"$bams1" | wc -l)
N_IND2=$(cat "$dir2"/"$bams2" | wc -l)
N_IND3=$(cat "$dir2"/"$bams3" | wc -l)
echo -e "Number of $taxa1 individuals: $N_IND1 \nNumber of $taxa2 individuals: $N_IND2 \nNumber of $taxa3 individuals: $N_IND3\n"
CHR1=$(( $N_IND1 * 2 ))
CHR2=$(( $N_IND2 * 2 ))
CHR3=$(( $N_IND3 * 2 ))
echo -e "Number of $taxa1 chromosomes: $CHR1 \nNumber of $taxa2 chromosomes: $CHR2 \nNumber of $taxa3 chromosomes: $CHR3\n"

echo "Calculating genotypes, posterior probabilities, SFS, etc..."
"$dir1"/angsd -bam "$dir2"/"$bams1" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND1" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa1"_r.indF -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa1"
"$dir1"/angsd -bam "$dir2"/"$bams2" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND2" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa2"_r.indF -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa2"
"$dir1"/angsd -bam "$dir2"/"$bams3" -anc "$dir4"/"$anc" -minMapQ "$minMapQ" -minQ "$minQ" -doMajorMinor 1 -doMaf 2 -doSaf 1 -GL 1 -minMaf 0.05 -P "$threads" -out "$taxa3".2pops
"$dir1"/angsd -bam "$dir2"/"$bams3" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND3" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa3"_r.indF -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa3".pca
echo "Done calculating genotypes, posterior probabilities, SFS files."
date

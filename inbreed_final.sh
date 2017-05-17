#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A-inbrd-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A-inbrd-stderr.txt
#SBATCH -J inbrd
#SBATCH -a 1
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 2
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

#Declare other variables
taxa="PA"						# file prefix
bams="angsd.admix.txt"				# file of BAM files for taxa
ref="Prunus_persica_v1.0_scaffolds.fa"			# file name for reference genome
anc="PC01.fa"						# file name for ancestral genome
threads="2"						# number of threads used
minMapQ="30"						# minimum mapping quality
minQ="20"						# mimimu base quality
#r="scaffold_""$x"":"					# region

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

# Utilized much of Simon's ANGSD workflow information on github
# https://github.com/XLEvolutionist/angsdworkflow/blob/master/first_teosinte_parents.html


##### CALCULATING GENOTYPE LIKELIHOODS FROM ANGSD #####
# Calculating the final genotype likelihoods because received following error message
# "Potential problem: -doMajorMinor 1 is based on genotype likelihoods,"
# "you must specify a genotype likelihood model -GL"
# Added -GL 1 like in genotype likelihood analyses without inbreeding
# Modified in full inbreeding script

date
echo "Calculating genotypes and posterior probabities with inbreeding..."
N_IND=$(cat "$dir2"/"$bams" | wc -l)
echo "Number of individuals: " $N_IND
N_SITES=$(( `zcat "$taxa"_F.mafs.gz | wc -l`-1 ))
echo "Number of sites: " $N_SITES
"$dir1"/angsd -bam "$dir2"/"$bams" -ref "$dir3"/"$ref" -anc "$dir4"/"$anc" -nInd "$N_IND" -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 2 -GL 1 -indF "$taxa"_F.indF -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa"_F_final
echo "Done calculating genotypes."
date

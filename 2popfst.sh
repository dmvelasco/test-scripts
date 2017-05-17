#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-2popfst-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-2popfst-stderr.txt
#SBATCH -J 2popfst
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 12
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
taxa3="PA"				# file prefix for combined dulcis and persica
bams1="angsd.dulcis2.txt"		# list of BAM files for dulcis
bams2="angsd.persica2.txt"		# list of BAM files for persica
bams3="angsd.domesticates.txt"		# list of BAM files for dulcis and persica combined
ref="Prunus_persica_v1.0_scaffolds.fa"	# file name for reference genome
anc="PC01.fa"				# file name for ancestral genome
threads="12"				# number of threads used
minMapQ="30"				# minimum mapping quality
minQ="20"				# mimimum base quality
r="scaffold"				# region prefix

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

##### SAF calculation
# First calculate per pop saf for each population
# already done previously, below is just an example
#date
#echo "Calculating site allele frequencies (SAFs)..."
#"$dir1"/angsd -b list1 -anc hg19ancNoChr.fa -out pop1 -dosaf 1 -gl 1
#"$dir1"/angsd -b list2 -anc hg19ancNoChr.fa -out pop2 -dosaf 1 -gl 1
#echo "Done calculating SAFs."

##### 2D SFS
# Using the new version of ANGSD (saf.idx files)
# In case of low coverage sequencing, an improvement in the estimation accuracy can be
# achieved by using the joint-SFS (2D-SFS) as a prior for the sample allele frequency
# posterior distributions, as shown in Fumagalli et al. (2013):

date
echo "Calculating the 2D site frequency spectrum (SFS)..."
# Calculate the 2dsfs prior
"$dir5"/realSFS "$taxa1"_F_final.saf.idx "$taxa2"_F_final.saf.idx -P "$threads" > "$taxa3".ml
date
echo "Done calculating the 2D SFS."

# If needed, this estimated 2D-SFS can be converted, using the script convert.2Dsfs.to.dadi.R

# Below from http://www.popgen.dk/angsd/index.php/SFS_Estimation
# Posterior of the per-site distributions of the sample allele frequency
# If you supply a prior for the SFS (which can be obtained from the -doSaf/realSFS
# analysis), the output of the .saf file will no longer be site allele frequency
# likelihoods but instead will be the log posterior probability of the sample allele
# frequency for each site in logspace.

# http://popgen.dk/angsd/index.php/2d_SFS_Estimation
# ANGSD wiki: https://github.com/ANGSD/angsd/wiki/SFS


##### FST
#http://www.popgen.dk/angsd/index.php/Fst_PCA (2015-06-17)
#New Fancy Version for two populations with real data
date
echo "Calculating Fst..."

# Prepare the fst for easy window analysis etc
"$dir5"/realSFS fst index "$taxa1"_F_final.saf.idx "$taxa2"_F_final.saf.idx -sfs "$taxa3".ml -fstout "$taxa3"

# Get the global FST estimate
# Produces the global unweighted and weighted FST
# found in the standard out file (maybe standard error file)
"$dir5"/realSFS fst stats "$taxa3".fst.idx

# According to ANGSD website: "below is not tested that much, but seems to work"
"$dir5"/realSFS fst stats2 "$taxa3".fst.idx -win 1000 -step 50 > "$taxa3"_windows

echo "Done calculating Fst."
date

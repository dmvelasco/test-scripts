#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-pca-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-pca-stderr.txt
#SBATCH -J pca
#SBATCH -p bigmemm
#SBATCH -a 1-8
#SBATCH -n 1
#SBATCH -c 1
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
taxa="PP"						# file prefix
bams="angsd.persicabam.txt"				# file of BAM files for taxa
ref="Prunus_persica_v1.0_scaffolds.fa"			# file name for reference genome
anc="PC01.fa"						# file name for ancestral genome
threads="4"						# number of threads used
minMapQ="30"						# minimum mapping quality
minQ="20"						# mimimu base quality
r="scaffold_""$x"":"                                    # region

# ngsPopGen sample scripts https://github.com/mfumagalli/ngsPopGen/tree/master/examples

date
echo "Calculating PCA"
N_IND=$(cat "$dir2"/"$bams" | wc -l)
N_SITES=$(( `zcat "$taxa"_s"$x".mafs.gz | wc -l`-1 ))
echo -e "Number of individual: $N_IND \nNumber of sites: $N_SITES"
date
echo "Calculating genotypes, posterior probabilities, SFS, etc..."
N_IND=$(cat "$dir2"/"$bams" | wc -l)
echo "$N_IND"
"$dir1"/angsd -bam "$dir2"/"$bams" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa"_r.indF -r "$r" -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa"_s"$x"
echo "Done calculating genotypes."
date
gunzip -c "$taxa"_s"$x".geno.gz > "$taxa"_s"$x".geno

# Principal Component Analysis (PCA)
# Covariance matrix

# Use ngsCovar to estimate a covariance matrix between pairs of individuals. Use -norm 0 for low coverage
# sequencing data, which disables normalization proposed in Patterson et al. (2006).
# The first way is to compute an approximation of the posterior of the covariance matrix, by weighting each
# site by its probability of being variable, as proposed in Fumagalli et al. (2013):
#../ngsCovar -probfile testA.geno -outfile testA.covar1 -nind 24 -nsites 10000 -call 0 -sfsfile testA.rf.saf -norm 0
#"$dir6"/ngsCovar -probfile "$taxa"_s"$x".geno.gz -outfile "$taxa"_s"$x".covar1 -nind "$N_IND" -nsites "$N_SITES" -call 0 -sfsfile "$taxa"_s"$x".saf -norm 0
# .saf.pos.gz or .saf?

# Alternatively, one can remove non-variable or low-frequency sites (this is the preferred way under most cases) with the option -minmaf:
#../ngsCovar -probfile testA.geno -outfile testA.covar2 -nind 24 -nsites 10000 -call 0 -minmaf 0.05
"$dir6"/ngsCovar -probfile "$taxa"_s"$x".geno -outfile "$taxa"_s"$x".covar -nind "$N_IND" -nsites "$N_SITES" -block_size 20000 -call 0 -minmaf 0.05

# Finally, in case of high sequencing depth, one can call genotypes as the genotype with the highest posterior probability:
#../ngsCovar -probfile testA.geno -outfile testA.covar3 -nind 24 -nsites 10000 -call 1 -minmaf 0.05
#"dir6"/ngsCovar -probfile testA.geno -outfile testA.covar3 -nind 24 -nsites 10000 -call 1 -minmaf 0.05

# These commands will produce text files with a symmetric covariance matrix NxN (for N individuals).

# PCA plot
# From the covariance matrix we can use a simple R script to perform an eigenvalue decomposition and plot the
# PCA results. Please make sure the use updated versions of required R packages. First, let's create a dummy plink cluster file.
#Rscript --vanilla --slave -e 'write.table(cbind(seq(1,24),rep(1,24),c(rep("A",10),rep("B",8),rep("C",6))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="testA.clst", quote=F)'

# Assuming we want to plot the first 2 PCA components from the testA.covar1 file, we can use the plotPCA.R script provided:
#Rscript --vanilla --slave $SCRIPTS/plotPCA.R -i testA.covar1 -c 1-2 -a testA.clst -o testA.pca.SAF.pdf

# Please note that you need 'ggplot2' and 'optparse' R libraries installed. This script will output the
# explained genetic variance for each component and save as output the PCA plot.
echo "Done calculating PCA"
date

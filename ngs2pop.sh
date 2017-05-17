#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-2pop-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-2pop-stderr.txt
#SBATCH -J 2pop
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
taxa3="PE"				# file prefix for combined dulcis and persica
bams1="angsd.dulcisbam.txt"		# list of BAM files for dulcis
bams2="angsd.persicabam.txt"		# list of BAM files for persica
bams3="angsd.domesticatesbam.txt"	# list of BAM files for dulcis and persica combined
ref="Prunus_persica_v1.0_scaffolds.fa"	# file name for reference genome
anc="PC01.fa"				# file name for ancestral genome
threads="8"				# number of threads used
minMapQ="30"				# minimum mapping quality
minQ="20"				# mimimum base quality
r="scaffold"				# region prefix

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

##### CALCULATING GENOTYPES, GENOTYPE LIKELIHOODS, SFS, MAFS, ETC. WHILE ACCOUNTING FOR INBREEDING

date
echo "Prepping variables..."
N_IND1=$(cat "$dir2"/"$bams1" | wc -l)
N_IND2=$(cat "$dir2"/"$bams2" | wc -l)
N_IND3=$(cat "$dir2"/"$bams3" | wc -l)
echo -e "Number of $taxa1 individuals: $N_IND1 \nNumber of $taxa2 individuals: $N_IND2 \nNumber of $taxa3 individuals: $N_IND3"
CHR1=$(( $N_IND1 * 2 ))
CHR2=$(( $N_IND2 * 2 ))
CHR3=$(( $N_IND3 * 2 ))
echo -e "Number of $taxa1 chromosomes: $CHR1 \nNumber of $taxa2 chromosomes: $CHR2 \nNumber of $taxa3 chromosomes: $CHR3"
date

##### 2D SFS
# New ANGSD 2dSFS method with realSFS
#date
#echo "Estimating 2D SFS..."

# Below from http://www.popgen.dk/angsd/index.php/SFS_Estimation
# Posterior of the per-site distributions of the sample allele frequency
# If you supply a prior for the SFS (which can be obtained from the -doSaf/realSFS
# analysis), the output of the .saf file will no longer be site allele frequency
# likelihoods but instead will be the log posterior probability of the sample allele
# frequency for each site in logspace.

# http://popgen.dk/angsd/index.php/2d_SFS_Estimation
# ANGSD wiki: https://github.com/ANGSD/angsd/wiki/SFS
# BELOW WORKS
#for i in {1..8}; do
#	"$dir5"/realSFS "$taxa1".saf.idx -P "$threads" -r "$r"_"$i": >> "$taxa1".saf.sfs
#	"$dir5"/realSFS "$taxa2".saf.idx -P "$threads" -r "$r"_"$i": >> "$taxa2".saf.sfs
#	"$dir5"/realSFS "$taxa1".saf.idx "$taxa2".saf.idx -P "$threads" -r "$r"_"$i": >> "$taxa3".2pops.2d.sfs
#	"$dir5"/realSFS "$taxa1".saf.idx -maxIter 100 -P "$threads"  -r "$r"_"$i": >> "$taxa1".sfs
#	"$dir5"/realSFS "$taxa2".saf.idx -maxIter 100 -P "$threads"  -r "$r"_"$i": >> "$taxa2".sfs
#done

# If needed, this estimated 2D-SFS can be converted, using the script convert.2Dsfs.to.dadi.R
#echo "Done estimating SFS for separate and combined populations
#date


##### FST
# Generate .saf files from each population using ANGSD SFS Estimation
#	using a 2D-SFS as a prior, estimated using ngs2dSFS
#	using marginal spectra as priors, estimated using realSFS
# tutorial at https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
# https://github.com/ANGSD/angsd/issues/14
###### If you are using the new version of angsd, the ones that generate
# saf.idx files, then you can simply run realSFS on the saf.idx files
# without having to rerun angsd to get the intersect.
# ./realSFS PD.saf.idx PP.saf.idx -P 24 > 2dsfs.txt
date
echo "Calculating Fst using 2dSFS as prior..."
Rscript /home/dmvelasc/Software/ngsTools/scripts/convertSFS.R "$taxa3".2pops.2d.sfs > "$taxa3".2d.sfs
# works now that path to R script added
N_SITES=$(cat "$taxa3".2d.sfs | wc -m)
echo -e "Number of sites: $N_SITES"
"$dir6"/ngsFST -postfiles "$taxa1".saf "$taxa2".saf -priorfile "$taxa3".2d.sfs -islog 1 -nind "$N_IND1" "$N_IND2" -nsites "$N_SITES" -outfile "$taxa3".fst
#added -islog

# In case of low coverage sequencing, an improvement in the estimation accuracy can be
# achieved by using the joint-SFS (2D-SFS) as a prior for the sample allele frequency
# posterior distributions, as shown in Fumagalli et al. (2013):
#echo "Done calculating Fst."

#######################################
# COMBINED DULCIS AND PERSICA SAMPLES #
#######################################

##### ADMIXTURE - WORKS
#date
#echo "Calculating admixture..."
#example for creating Beagle input file, -doGlf 2; use combined bam file list
#"$dir1"/angsd -GL 1 -out "$taxa3".admix -nThreads "$threads" -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam "$dir2"/"$bams3"
#"$dir1"/NGSadmix -likes "$taxa3".admix.beagle.gz -K 2 -P "$threads" -o "$taxa3"_admix -minMaf 0.05
#echo "Done calculating admixture."
#date
########

##### PRINCIPAL COMPONENT ANALYSIS (PCA) - testing separately
# Covariance matrix
#date
#echo "Calculating PCA"
#gunzip -c "$taxa3".pca.geno.gz > "$taxa3".pca.geno

# Use ngsCovar to estimate a covariance matrix between pairs of individuals. Use -norm 0
# for low coverage sequencing data, which disables normalization proposed in Patterson et
# al. (2006). The first way is to compute an approximation of the posterior of the
# covariance matrix, by weighting each site by its probability of being variable, as
# proposed in Fumagalli et al. (2013):
#../ngsCovar -probfile testA.geno -outfile testA.covar1 -nind 24 -nsites 10000 -call 0 -sfsfile testA.rf.saf -norm 0

# Alternatively, one can remove non-variable or low-frequency sites (this is the preferred way under most cases) with the option -minmaf:
#../ngsCovar -probfile testA.geno -outfile testA.covar2 -nind 24 -nsites 10000 -call 0 -minmaf 0.05
#N_SITES2=$( `zcat "$taxa3".pca.mafs.gz | wc -l`-1 )
#"$dir6"/ngsCovar -probfile "$taxa3".pca.geno -outfile "$taxa3".covar -nind "$N_IND3" -nsites "$N_SITES2" -block_size 20000 -call 0 -minmaf 0.05

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
#echo "Done calculating PCA"
#date

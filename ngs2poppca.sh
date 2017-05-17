#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-pca-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-pca-stderr.txt
#SBATCH -J 2poppca
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 4
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
dir6="/home/dmvelasc/Software/ngsPopGen"		# ngsPopGen directory #formerly used in subdirectory of ngsTools

#Declare other variables
taxa1="PD"				# file prefix for dulcis
taxa2="PP"				# file prefix for persica
taxa3="PE"				# file prefix for combined dulcis and persica
#taxa3="PA"				# file prefix for combined dulcis and persica
bams1="angsd.dulcis2.txt"		# list of BAM files for dulcis
bams2="angsd.persica2.txt"		# list of BAM files for persica
bams3="angsd.domesticates.txt"		# list of BAM files for dulcis and persica combined
#bams3="angsd.admix.txt"		# list of BAM files for dulcis and persica combined
ref="Prunus_persica_v1.0_scaffolds.fa"	# file name for reference genome
anc="PC01.fa"				# file name for ancestral genome
threads="4"				# number of threads used
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


#######################################
# COMBINED DULCIS AND PERSICA SAMPLES #
#######################################

##### PRINCIPAL COMPONENT ANALYSIS (PCA)
# Covariance matrix
date
echo "Calculating PCA data"
#BELOW works
"$dir1"/angsd -bam "$dir2"/"$bams3" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND3" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa3"_F.indF -r scaffold_1 -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa3"_s1.pca
# without region
#"$dir1"/angsd -bam "$dir2"/"$bams3" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND3" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa3"_F.indF -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa3".pca
# entire genome data was too big to unzip without eating all drive space so limiting to chr 1
#BELOW works, once done not needed again. should probably have a whether file exists check.
echo "Unzipping geno.gz file."
gunzip -c "$taxa3"_s1.pca.geno.gz > "$taxa3"_s1.pca.geno
#gunzip -c "$taxa3".pca.geno.gz > "$taxa3".pca.geno

# Use ngsCovar to estimate a covariance matrix between pairs of individuals. Use -norm 0
# for low coverage sequencing data, which disables normalization proposed in Patterson et
# al. (2006). The first way is to compute an approximation of the posterior of the
# covariance matrix, by weighting each site by its probability of being variable, as
# proposed in Fumagalli et al. (2013):
#../ngsCovar -probfile testA.geno -outfile testA.covar1 -nind 24 -nsites 10000 -call 0 -sfsfile testA.rf.saf -norm 0
N_SITES2=$(( `zcat "$taxa3"_s1.pca.mafs.gz | wc -l`-1 ))
#N_SITES2=$(( `zcat "$taxa3".pca.mafs.gz | wc -l`-1 ))
echo $N_SITES2
# below test based on responses to issue on ngsPopGen github repo
"$dir5"/realSFS print "$taxa3"_s1.pca.saf.idx -oldout 1
#"$dir5"/realSFS print "$taxa3".pca.saf.idx -oldout 1

#use below line while testing otherwise ngsCovar errors that there is already a file by that name
#rm "$taxa3"_s1_wtd.covar
echo "Calculating weighted covariance # not this time"
#"$dir6"/ngsCovar -probfile "$taxa3"_s1.pca.geno -outfile "$taxa3"_s1_wtd.covar -nind "$N_IND3" -nsites "$N_SITES2" -block_size 20000 -call 0 -sfsfile "$taxa3"_s1.pca.saf -norm 0
#"$dir6"/ngsCovar -probfile "$taxa3".pca.geno.gz -outfile "$taxa3"_wtd.covar -nind "$N_IND3" -nsites "$N_SITES2" -block_size 20000 -call 0 -sfsfile "$taxa3".pca.saf -norm 0

# -sfsfile	"file with sample allele frequency posterior probabilities"
#		The example in the main README shows a file named "pop.sfs.ml.norm" for this parameter (2014-01-03)
#		The examples directory example shows a file named "testA.rf.saf" for this parameter (2014-01-07)
#		The error when using PE.ml generated for another analysis is "Possible error, binary files might be broken",
#		which suggests the testA.rf.saf file may be the correct one
# -norm		"0 [no normalization, recommended if no SNP calling is performed]"
#		"1 [matrix is normalized by p(1-p) as in Patterson et al. (2006)]"
#		"2 [normalized by 2p(1-p)]"

# NOTE: 2015-05-19 commit "Added support to read GENO and SFS files directly from STDIN and print results directly to STDOUT"

# Alternatively, one can remove non-variable or low-frequency sites (this is the preferred way under most cases) with the option -minmaf:
#../ngsCovar -probfile testA.geno -outfile testA.covar2 -nind 24 -nsites 10000 -call 0 -minmaf 0.05

#DV NOTE: works, but gives odd results for samples with low coverage.
echo "Calculating covariance with min MAF of 0.05"
#N_SITES2=$(( `zcat "$taxa3"_s1.pca.mafs.gz | wc -l`-1 ))
"$dir6"/ngsCovar -probfile "$taxa3"_s1.pca.geno -outfile "$taxa3"_s1.covar -nind "$N_IND3" -nsites "$N_SITES2" -block_size 20000 -call 0 -minmaf 0.05
#"$dir6"/ngsCovar -probfile "$taxa3".pca.geno.gz -outfile "$taxa3".covar -nind "$N_IND3" -nsites "$N_SITES2" -block_size 20000 -call 0 -minmaf 0.05


# Finally, in case of high sequencing depth, one can call genotypes as the genotype with the highest posterior probability:
#../ngsCovar -probfile testA.geno -outfile testA.covar3 -nind 24 -nsites 10000 -call 1 -minmaf 0.05
#"dir6"/ngsCovar -probfile testA.geno -outfile testA.covar3 -nind 24 -nsites 10000 -call 1 -minmaf 0.05

# These commands will produce text files with a symmetric covariance matrix NxN (for N individuals).
echo "Done calculating PCA"
date

# PCA plot
# From the covariance matrix we can use a simple R script to perform an eigenvalue decomposition and plot the
# PCA results. Please make sure the use updated versions of required R packages. First, let's create a dummy plink cluster file.
#Rscript --vanilla --slave -e 'write.table(cbind(seq(1,24),rep(1,24),c(rep("A",10),rep("B",8),rep("C",6))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="testA.clst", quote=F)'

# Assuming we want to plot the first 2 PCA components from the testA.covar1 file, we can use the plotPCA.R script provided:
#Rscript --vanilla --slave $SCRIPTS/plotPCA.R -i testA.covar1 -c 1-2 -a testA.clst -o testA.pca.SAF.pdf

# Please note that you need 'ggplot2' and 'optparse' R libraries installed. This script will output the
# explained genetic variance for each component and save as output the PCA plot.

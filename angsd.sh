#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-PPtajd-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-PPtajd-stderr.txt
#SBATCH -J angsd
#SBATCH -p bigmemm
#SBATCH -a 1-8
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
taxa="PD"						# file prefix
bams="angsd.dulcis2.txt"			# file of BAM files for taxa
ref="Prunus_persica_v1.0_scaffolds.fa"			# file name for reference genome
anc="PC01.fa"						# file name for ancestral genome
threads="4"						# number of threads used
minMapQ="30"						# minimum mapping quality
minQ="20"						# mimimu base quality
r="scaffold_""$x"":"                                    # region

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

# Utilized much of Simon's ANGSD workflow information on github
# https://github.com/XLEvolutionist/angsdworkflow/blob/master/first_teosinte_parents.html

N_IND=$(cat "$dir2"/"$bams" | wc -l)
echo "$N_IND"
CHROM=$(( N_IND*2 ))
echo "$CHROM"

# Based on popgen.dk/angsd/index.php/Tajima

##### CALCULATING GENOTYPES, GENOTYPE LIKELIHOODS, SFS, MAFS, ETC. WHILE ACCOUNTING FOR INBREEDING
date
echo "Calculating genotypes, posterior probabilities, SFS, etc..."
"$dir1"/angsd -bam "$dir2"/"$bams" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -nInd "$N_IND" -doMajorMinor 1 -doMaf 1 -doSaf 2 -GL 1 -doPost 1 -doGeno 32 -indF "$taxa"_F.indF -r "$r" -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa"_"$x"
echo "Done calculating genotypes."
date

# -bam		"$dir2"/"$bams"			Path to list of BAM files for analysis
# -doSaf	2				site allele frequency,
# -out		"$taxa"				Prefix for outfiles
# -anc		"$dir4"/"$anc"			Path of ancestral genome sequence file
# -ref		"$dir3"/"$ref"			Path of reference genome sequence file
# -GL		1				Genotype likelihood
# -nThreads	"$threads"			number of threads
# -rf		"$dir2"/angsd.regions.txt	Path of file containing regions for analysis
# -indF		"$taxa".approx_indF		Inbreeding coefficient
# -doMaf	1				x
# -doMajorMinor	1				x
# -minMapQ	30				mimimum mapping quality
# -minQ		20				minimum base quality


##### CALCULATING THE SITE FREQUENCY SPECTRUM
# The .saf file generated above is needed first.
# repeats all of above, except GL so added above
# You need to declare the number of chromosomes (so for 10 individuals, you choose 20, and for 20 you choose 40).

## this is done below
echo "Calculating the SFS: "
date
# 1D SFS (NEW - WORKS, errors if number of chromosomes is included, probably due to the newer ANGSD version)
"$dir5"/realSFS "$taxa"_"$x".saf.idx -P "$threads" -r "$r" -maxIter 1000 > "$taxa"_"$x".sfs
echo "Done calculating SFS "
date

# "$taxa".saf			per site allele frequency file from SFS estimation
# $CHROM	"$CHROM"	haploid number of sample chromosomes
# -maxIter	100		maximum number of iterations for EM algorithm
# -P		"$threads"	number of threads
# "$taxa".sfs			SFS output file


##### ESTIMATING GENOME- OR CHROMOSOME-WIDE THETAS #####
# Calculate the genome wide thetas with the following:

date
echo "Estimating thetas... "
"$dir1"/angsd -bam "$dir2"/"$bams" -anc "$dir4"/"$anc" -ref "$dir3"/"$ref" -doThetas 1 -doSaf 1 -pest "$taxa"_"$x".sfs -GL 2 -P "$threads" -r "$r" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa"_"$x".sfs
echo "Done estimating thetas. "
date

# -bam		"$dir2"/"$bams"			bam file paths
# -out		"$taxa".sfs			out file
# -doThetas	1				x
# -doSaf	1				x
# -pest		"$taxa".sfs			p estimate
# -anc		"$dir4"/"$anc"			ancestral genome
# -ref		"$dir3"/"$ref"			reference genome
# -GL		2				genotype likelihood 2
# -P		"$threads"			number of threads
# -r		"$r"				region to analyze
# -minMapQ	"$minMapQ"			minimum mapping quality
# -minQ		"$minQ"				minimum base quality


##### CALCULATING TAJIMA'S D #####
# Remove "scaffold" from sfs.thetas.gz file.
# Make a BED file of the genome-wide theta values and calculate Tajima’s D.

date
echo "Calculating Tajima's D... "
# perl to remove "scaffold" from file so chromosome is just number
zcat "$taxa"_"$x".sfs.thetas.gz | perl -plne 's/scaffold_(\w+)/$1/' - > "$taxa"_"$x"_mod.sfs.thetas.gz
# create a binary version of the .thetas.gz
"$dir5"/thetaStat make_bed "$taxa"_"$x"_mod.sfs.thetas.gz
# calculate Tajimas D
"$dir5"/thetaStat do_stat "$taxa"_"$x"_mod.sfs.thetas.gz -nChr "$CHROM" -r "$x" -win 1000 -step 50 -outnames "$taxa"_"$x".thetasW1kS50
echo "Done calculating Tajima's D."
date

# Line 2: do_stat
# out.thetas.gz
# -nChr		16
# -win		1000
# -step		50
# -outnames	"$taxa"_x_thetasW1kS50


####### CALCULATING SUMMARY STATISTICS
date
echo "Calculating Summary Statistics..."
N_SITES=$(( `zcat "$taxa"_"$x"_tajd.mafs.gz | wc -l`-1 ))
echo "$N_SITES"
"$dir6"/ngsStat -npop 1 -postfiles "$taxa"_"$x".saf -nsites "$N_SITES" -iswin 1 -nind "$N_IND" -outfile "$taxa"_"$x".stat -isfold 0 -islog 1 -block_size "$N_SITES"

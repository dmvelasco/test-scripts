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
# Calculate genotype likelihoods for each sample in order to estimate the inbreeding coefficient (F)

date
echo "Calculating genotype likelihoods..."
N_IND=$(cat "$dir2"/"$bams" | wc -l)
echo "Number of individuals: " $N_IND
"$dir1"/angsd -bam "$dir2"/"$bams" -ref "$dir3"/"$ref" -anc "$dir4"/"$anc" -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -GL 1 -doGLF 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-4 -out "$taxa"_F
echo "Done calculating genotype likelihoods."
date

# -bam		<file>		list of bam file names and paths
# -doGLF	3		flag to output binary format genotype likelihood file
#				this is the file format needed for ngsF to estimate F
# -GL		1		calculate genotype likelihoods using SAMtools algorithm
#				details at "https://www.broadinstitute.org/gatk/media/docs/Samtools.pdf"
# -out		"$taxa"		prefix for outfiles
# -doMaf	2		Calculate per site frequencies; 2 is fixed major unknown minor; 1 is fixed major and minor
# -SNP_pval	1e-6		For polymorphic sites only work with sites with a p-value less than [float], requires -doMaf
# -doMajorMinor	1		Infer the major and minor alleles, needed to estimate allele frequencies from genotype likelihoods
# -nThreads	"$threads"	number of threads used for computation
# -r		"$r"		region for analysis
# -minMapQ	30		minimum mapping quality
# -minQ		20		minimum base quality

# From https://github.com/fgvieira/ngsF/tree/master/examples
# In this example, we will estimate inbreeding coefficients per individual and incorporate them into
# the calculation of posterior probabilities. First, calculate genotype likelihoods and call SNPs:
#$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd 20 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-4 -out testF.HWE


##### INBREEDING COEFFICIENTS FROM NGSF ######
# Estimating the inbreeding coefficient (F) for each sample

# The following produces two output files. The first line estimates F and then these estimates
# are parsed to the second line which then refines the estimates. These estimates can then be
# combined with their samples names using a unix command like the following:
# 	paste "$taxa"_F.approx_indF "$dir2"/"$bams" > "$taxa"_F.approx_named_indF.txt
# where "$bams" is the list of .bam files used in the analysis.

date
echo "Estimating inbreeding coefficients..."
gunzip -c "$taxa"_F.glf.gz > "$taxa"_F.glf
N_SITES=$(( `zcat "$taxa"_F.mafs.gz | wc -l`-1 ))
# double parenthesis indicate mathematic operation
# backticks indicate executing the command first
# the lines in the file is read, lines are counted then reduced by 1 b/c of header line in mafs.gz
echo "Number of sites: " $N_SITES
"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_F.glf --out "$taxa"_F.approx_indF --approx_EM --seed 12345 --init_values r --n_threads "$threads"
"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_F.glf --out "$taxa"_F.indF --init_values "$taxa"_F.approx_indF.pars --n_threads "$threads"

# PREVIOUS TESTS, above that are _F were _r
#"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-3 --glf "$taxa".glf --out "$taxa"_e.approx_indF --approx_EM --seed 12345 --init_values e --n_threads "$threads"
#"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-3 --glf "$taxa".glf --out "$taxa"_e.indF --init_values "$taxa"_e.approx_indF.pars --n_threads "$threads"
#"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-3 --glf "$taxa".glf --out "$taxa"_ee.indF --init_values e --n_threads "$threads"

#"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-3 --glf "$taxa".glf --out "$taxa"_rr.indF --init_values r --n_threads "$threads"
# added -- in front of functions due to recently updated examples at ngsF github repo https://github.com/fgvieira/ngsF
echo "Done estimating inbreeding coefficients."
date

# Common to both lines
# -n_ind	$N_IND				Sample size (number of individuals), variable created from
#						reading the list of used BAMS
# -n_sites	$N_SITES			Total number of sites, variable created from reading MAF file
# -min_epsilon 	0.001				Maximum RMSD between iterations to assume convergence [1e-5]
# -glf 		"$taxa".glf			Input uncompressed genotype likelihood file
# -n_threads	"$threads"			Number of threads used for analysis [1]
# -out		"$taxa".approx_indF		Output file name
# -init_values	line 1 = r			Initiatl values of individual F and site frequency. Can be
#		line 2 = file			(r)andom, (e)stimated from data, (u)niform at 0.01, or read
#						from a file

# Specific to first line
# -approx_EM					Use the faster approximated EM ML algorithm
# -seed		12345				not defined at github repo

date
echo "Calculating genotypes and posterior probabilities..."
"$dir1"/angsd -bam "$dir2"/"$bams" -ref "$dir3"/"$ref" -anc "$dir4"/"$anc" -nInd "$N_IND" -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 2 -GL 1 -indF "$taxa"_F.indF -nThreads "$threads" -minMapQ "$minMapQ" -minQ "$minQ" -out "$taxa"_F_final
echo "Done calculating genotypes."
date

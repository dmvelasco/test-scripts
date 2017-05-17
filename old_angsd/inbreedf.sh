#!/bin/bash -l
#OUTDIR=/home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/angsd_output
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-inbrd-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-inbrd-stderr.txt
#SBATCH -J inbrd
#SBATCH -a 1-8
#SBATCH -p bigmemm
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

#Declare other variables
taxa="PD"						# file prefix
bams="angsd.dulcisbam.txt"				# file of BAM files for taxa
ref="Prunus_persica_v1.0_scaffolds.fa"			# file name for reference genome
anc="PC01.fa"						# file name for ancestral genome
threads="1"						# number of threads used
minMapQ="30"						# minimum mapping quality
minQ="20"						# mimimu base quality
r="scaffold_""$x"":"					# region

# ANGSD resource at http://cgrlucb.wikispaces.com/file/view/CGRL_SNP_workshop.pdf

# Utilized much of Simon's ANGSD workflow information on github
# https://github.com/XLEvolutionist/angsdworkflow/blob/master/first_teosinte_parents.html


##### CALCULATING GENOTYPE LIKELIHOODS FROM ANGSD #####
# Calculate genotype likelihoods for each sample in order to estimate the inbreeding coefficient (F)

date
echo "Calculating genotype likelihoods..."
N_IND=$(cat "$dir2"/"$bams" | wc -l)
echo "Number of individuals: " $N_IND
"$dir1"/angsd -bam "$dir2"/"$bams" -doGlf 3 -GL 1 -out "$taxa"_s"$x" -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doGeno 9 -doPost 1 -ref "$dir3"/"$ref" -anc "$dir4"/"$anc" -nThreads "$threads" -r "$r" -minMapQ "$minMapQ" -minQ "$minQ"
echo "Done calculating genotype likelihoods."
date

# -bam		<file>		list of bam file names and paths
# -doGLF	3		flag to output binary format genotype likelihood file
#				this is the file format needed for ngsF to estimate F
# -GL		1		calculate genotype likelihoods using SAMtools algorithm
#				details at "https://www.broadinstitute.org/gatk/media/docs/Samtools.pdf"
# -out		"$taxa"_s"$x"	prefix for outfiles with scaffold number
# -doMaf	2		Calculate per site frequencies; 2 is fixed major unknown minor; 1 is fixed major and mionr
# -SNP_pval	1e-6		For polymorphic sites only work with sites with a p-value less than [float], requires -doMaf
# -doMajorMinor	1		Infer the major and minor alleles, needed to estimate allele frequencies from genotype likelihoods
# -nThreads	"$threads"	number of threads used for computation
# -r		"$r"		region for analysis
# -minMapQ	30		minimum mapping quality
# -minQ		20		minimum base quality


##### INBREEDING COEFFICIENTS FROM NGSF ######
# Estimating the inbreeding coefficient (F) for each sample

# The following produces two output files. The first line estimates F and then these estimates
# are parsed to the second line which then refines the estimates. These estimates can then be
# combined with their samples names using a unix command like the following:
# 	paste "$taxa".approx_indF "$dir2"/"$bams" > "$taxa".approx_named_indF.txt
# where "$bams" is the list of .bam files used in the analysis.

echo ""

date
echo "Estimating inbreeding coefficients..."
gunzip -c "$taxa"_s"$x".glf.gz > "$taxa"_s"$x".glf
N_SITES=$(( `zcat "$taxa"_s"$x".mafs.gz | wc -l`-1 ))
# double parenthesis indicate mathematic operation
# backticks indicate executing the command first
# the lines in the file is read, lines are counted then reduced by 1 b/c of header line in mafs.gz
echo "Number of sites: " $N_SITES
"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_s"$x".glf --out "$taxa"_s"$x"_r.approx_indF --approx_EM --seed 12345 --init_values r --n_threads "$threads"
"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_s"$x".glf --out "$taxa"_s"$x"_r.indF --init_values "$taxa"_s"$x"_r.approx_indF.pars --n_threads "$threads"

"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_s"$x".glf --out "$taxa"_s"$x"_e.approx_indF --approx_EM --seed 12345 --init_values e --n_threads "$threads"
"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_s"$x".glf --out "$taxa"_s"$x"_e.indF --init_values "$taxa"_s"$x"_e.approx_indF.pars --n_threads "$threads"
"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_s"$x".glf --out "$taxa"_s"$x"_ee.indF --init_values e --n_threads "$threads"

"$dir1"/ngsF --n_ind $N_IND --n_sites $N_SITES --min_epsilon 1e-6 --glf "$taxa"_s"$x".glf --out "$taxa"_s"$x"_rr.indF --init_values r --n_threads "$threads"
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


##### ESTIMATING THE SITE FREQUENCY SPECTRUM #####
# Begin estimating the site frequency spectrum. The .saf file generated below is needed first.

#date
#echo "Estimating the SFS..."
#"$dir1"/angsd -bam "$dir2"/"$bams" -doSaf 2 -out "$taxa"_s"$x" -ref "$dir3"/"$ref" -anc "$dir4"/"$anc" -GL 1 -P 12 -r "$r" -indF "$taxa"_s"$x".approx_indF -doMaf 1 -doMajorMinor 1 -minMapQ "$minMapQ" -minQ "$minQ"
#echo "Done estimating the SFS."
#date

# -bam		"$dir2"/"$bams"			Path to list of BAM files for analysis
# -doSaf	2				
# -out		"$taxa"				Prefix for outfiles
# -anc		"$dir4"/"$anc"			Path of ancestral genome sequence file
# -ref		"$dir3"/"$ref"			Path of reference genome sequence file
# -GL		1				Genotype likelihood
# -P		12				
# -r		"$r"				Region for analysis
# -indF		"$taxa".approx_indF		Inbreeding coefficient
# -doMaf	1				
# -doMajorMinor	1				
# -minMapQ	30				mimimum mapping quality
# -minQ		20				minimum base quality


# You need to declare the number of chromosomes (so for 10 individuals, you choose 20, and for 20 you choose 40).

#date
#echo "Calculating SFS..."
#CHROM=$(( N_IND*2 ))
#echo $CHROM
#"$dir5"/realSFS "$taxa".saf $CHROM -maxIter 100 -P 4 > "$taxa".sfs
#echo "Done calculating SFS."
#date

# "$taxa".saf	per site allele frequency file from SFS estimation
# $CHROM	haploid number of sample chromosomes
# -maxIter	100	maximum number of iterations for EM algorithm
# -P		4	number of threads
# "$taxa".sfs	SFS output file

#####################
# BELOW ARE NOT QUITE WORKING, SWITCH TO USING ANGSD-WRAPPER
####################

##### ESTIMATING GENOME- OR CHROMOSOME-WIDE THETAS #####
# Calculate the genome wide thetas with the following:

#date
#echo "Estimating thetas..."
#"$dir1"/angsd -bam "$dir2"/"$bams" -out "$taxa".sfs -doThetas 1 -doSaf 1 -pest "$taxa".sfs -anc "$dir4"/"$anc" -GL 2 -P 4 -r "$r" -minMapQ "$minMapQ" -minQ "$minQ"
#echo "Done estimating thetas."
#date

# -bam		"$dir2"/"$bams"			
# -out		"$taxa".sfs			
# -doThetas	1				
# -doSaf	1				
# -pest		"$taxa".sfs			
# -anc		"$dir4"/"$anc"			
# -GL		2				
# -P		12				
# -r		"$r"				
# -minMapQ	"$minMapQ"			
# -minQ		"$minQ"				


##### CALCULATING TAJIMA'S D #####
# Make a BED file of the genome-wide theta values and calculate Tajima’s D:

# check if file extensions are correct

#date
#echo "Calculating Tajima's D..."
# create a binary version of thete.thetas.gz
#"$dir5"/thetaStat make_bed "$taxa".sfs.thetas.gz
# calculate Tajimas D
#"$dir5"/thetaStat do_stat "$taxa".sfs.thetas.gz -r "$r" -nChr 16 -win 1000 -step 200 -outnames "$taxa".thetasWin1KStep200_chr8.gz
#echo "Done calculating Tajima's D."
#date

# Line 2: do_stat
# out.thetas.gz
# -nChr		16
# -win		1000
# -step		200
# -outnames	"$taxa"thetasWindow_chr8.gz

##### Estimate coverage per locus per individual #####
# The purpose of this aspect of the analysis is to generate estimates of coverage
# over specific base pairs. This will allow us to assess the effects on SFS, pi
# and Tajima’s D depending on the number of samples that are available.
#echo "Starting Job: "
#date
#"$dir1"/angsd -bam "$dir2"/"$bams" -doCounts 1 -minInd 0 -dumpCounts 2 -minQ 20 -P 12 -r 10 -out "$taxa"_Chr10counts
#echo "Ending Job: "
#date


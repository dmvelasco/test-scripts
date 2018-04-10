#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemh
#SBATCH -t 4:00:00
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=64G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft)
module load conda3
#module load vcftools
#module load tabix

#############################
### Set up the parameters ###
#############################

####### ARRAYS #######

##############################
##### ACQUIRE SETUP DATA #####
## TWO FOR SPLIT ##
# path to sample list
list="/home/dmvelasc/Projects/Prunus/Script/smcpp_data.txt"
# view list to select pops/subpops, important below
# pop	subset	region			line number(s)
# PD	sub1-10				lines 1-10
# PD	sub1	range			line 1
# PD	sub7	China			line 7
# PD	sub8	C Asia			line 8
# PD	sub10	S Europe/W Asia		line 10
# PP	sub1-7				lines 11-17
# PP	sub1	range			line 11
# PP	sub3	China			line 13
# PP	sub4	China/Korea		line 14
# PM	sub1-3				lines 18-20
# PV	sub1-2				lines 21-22
# PS	all				line 23
# PG	all				line 24

line1=1		#line number of first pop/subpop wanted
line2=21	#line number of second pop/subpop wanted

# mapfile to extract first pop/subpop information and send to an array
mapfile -s "$line1" -n 1 -t id < "${list}"
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop1="${arr[0]}"
sub1="${arr[1]}"
samples1="${arr[2]}"

# mapfile to extract second pop/subpop information and send to an array
mapfile -s "$line2" -n 1 -t id < "${list}"
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop2="${arr[0]}"
sub2="${arr[1]}"
samples2="${arr[2]}"
################################

####### PATHS #######
# genome reference file
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# Joint VCF file - test joint VCF file
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_jointcalls.vcf"
# filtered joint VCF file - dulcis test VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_dulcis.recode.vcf"
# SMC++ prepped file
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"
#final="mu1_38x-8"
# smc file in $vcf_filt directory
smc_file="smcpp_prunus_biallelic.recode.vcf.gz"


####### PARAMETERS #######
mu="7.77e-9"	# population mutation rate
cut="5000"		# cutoff length for homozygosity

####################
### Begin script ###
####################

mkdir -p smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"

################## SPLIT #####################
echo -e "SMC++: create reciprocal joint frequency spectrum datasets\n and run split analysis"
date

# Create reciprocal joint frequency spectra
for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$smc_file" \
smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/"$pop1"_"$sub1"-"$pop2"_"$sub2"-"$i".smc.gz scaffold_"$i" \
"$pop1":"$samples1" "$pop2":"$samples2"
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$smc_file" \
smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/"$pop2"_"$sub2"-"$pop1"_"$sub1"-"$i".smc.gz scaffold_"$i" \
"$pop2":"$samples2" "$pop1":"$samples1"
done

# Refine the marginal estimates
smc++ split -o split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/ \
 smc_analysis/"$pop1"_"$sub1"_"$mu"/model.final.json \
 smc_analysis/"$pop2"_"$sub2"_"$mu"/model.final.json \
 smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"\*.smc.gz

##### posterior #####
# export (and visualize) the posterior distributon of the TMRCA

#smc++ posterior [--unfold | --polarization-error p] [--start START] [--end END]
# [--thinning k] [--heatmap heatmap.(pdf|png|gif|jpeg)] [--colorbar] [--M M]
# model.final.json arrays.npz data.smc[.gz] [data.smc[.gz] ...]

#positional arguments:
#  model.final.json      SMC++ model to use in forward-backward algorithm
#  arrays.npz            location to save posterior decoding arrays
#  data.smc[.gz]         SMC++ data set(s) to decode

#optional arguments:
#  -h, --help            show this help message and exit
#  -v, --verbose         increase debugging output, specify multiply times for more
#  --unfold              use unfolded SFS (alias for -p 0.0)
#  --polarization-error p, -p p
#                        uncertainty parameter for polarized SFS: observation (a,b) has probability [(1-p)*CSFS_{a,b} +
#                        p*CSFS_{2-a,n-2-b}]. default: 0.5
#  --start START         base at which to begin posterior decode
#  --end END             base at which to end posterior decode
#  --thinning k          emit full SFS only every <k>th site. default: 1
#  --heatmap heatmap.(pdf|png|gif|jpeg)
#                        Also draw a heatmap of the posterior TMRCA.
#  --colorbar            If plotting, add a colorbar

#HMM parameters:
#  --M M                 number of hidden states

################## SPLIT #####################

echo -e "SMC++: plot the joint analyses"
date

# Plot by generations
smc++ plot -c \
smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu".pdf \
smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/model.final.json \
smc_analysis/"$pop1"_"$sub1"_"$mu"/model.final.json \
smc_analysis/"$pop2"_"$sub2"_"$mu"/model.final.json

# Plot by years
smc++ plot -g 10 \
smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"_years.pdf \
smc_analysis/split/"$pop1"-"$line1"_"$pop2"-"$line2"_"$mu"/model.final.json \
smc_analysis/"$pop1"_"$sub1"_"$mu"/model.final.json \
smc_analysis/"$pop2"_"$sub2"_"$mu"/model.final.json

#-g	sets generation time in years to scale x-axis, otherwise in coalescent units
#--logy	plots the y-axis on a log scale
#-c	produces CSV-formatted table containing the data used to generate the plot

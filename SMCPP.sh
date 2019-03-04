#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemh
#SBATCH -t 24:00:00
#SBATCH -a 23-33%2
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=80G
set -e
set -u

####################
#   M O D U L E S  #
####################
# Add personal madule directory
module use /home/dmvelasc/MyModules

# Load personal and system modules
module load miniconda3
module load vcftools
module load tabix

###########################
#  E N V I R O N M E N T  #
###########################
# Activate python environment with SMC++
# set +u and -u is a work around for an unbound variable error in the activate script
set +u
source /home/dmvelasc/.virtualenvs/smcpp/bin/activate
set -u

#######################
#  V A R I A B L E S  #
#######################
## ARRAY VARIABLES ##
# Create array variables needed in script
x=$SLURM_ARRAY_TASK_ID
g=$(( x-1 ))

## PATHS AND FILES ##
# Create path variables used in script
# genome reference file
#genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# Joint VCF file
#vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/all_jointcalls.vcf"
# filtered joint VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# smc file in $vcf_filt directory
smc_file="smcpp_prunus_biallelic.recode.vcf.gz"
# SMC++ prepped file directory
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"
# sample list
list="/home/dmvelasc/Projects/Prunus/Script/smcpp_data.txt"


## ANALYSIS PARAMETERS ##
# estimate type [estimate (est) or cross validation (cv)]
type="est"
# population mutation rate, below are rates calculated in Xie et al. 2016
mu="7.77e-9"
  # Xie et al. 2016
  # 7.77e-9 (parent to selfed progeny)
  # 9.48e-9 (low heterozygosity peach to progeny)
  # 1.38e-8 (high heterozygosity peach to interspecific cross to selfed progeny)

# use one of the below for quality control for better estimation
# cutoff length for homozygosity, ignores runs greater than this length
#cut="5000"
# BED format mask file
mask_file="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.softmasked.bed.bgz"

## POPULATION, SUBPOPULATION, AND SAMPLE INFORMATION #####
# mapfile to extract sample ID and read name information, each line is array item
mapfile -s "$g" -n 1 -t id < "${list}"
# -s number of rows to skip; -n number of rows to read; -t (remove leading/trailing whitespace?)
# id is the array name (anything in this position is the array name)

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop="${arr[0]}"
sub="${arr[1]}"
samples="${arr[2]}"


##########################################################
#  B E G I N N I N G  O F  S C R I P T  C O M M A N D S  #
##########################################################

#echo -e "begin file preparation for SMC++ run for loop by chromosome:\n select individuals and populations at this step"
#date

##### SMC++ FINAL PREP #####
# only really need to do this part once
# convert files from VCF to SMC format
# *.smc.gz files for each chromosome are the SMC++ output file
#for i in {1..8}; do
# first runs used a cutoff value, found in directory ${smc_in}_1
#  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$smc_file" "$smc_in"/"${pop}_${sub}-${i}".smc.gz scaffold_"$i" "${pop}:${samples}"
# rerunning with v 1.15 used mask file
#  smc++ vcf2smc --mask "$mask_file" "$vcf_filt"/"$smc_file" "$smc_in"/"${pop}_${sub}-${i}".smc.gz scaffold_"$i" "${pop}:${samples}"
#done

mkdir -p smc_analysis/"${mu}"/"$type"/"${pop}_${sub}"

echo -e "begin SMC++ analysis"
date

##### SMC++ CROSS VALIDATION ANALYSIS #####
# model.final.json and iterations are outputs
smc++ estimate --cores 12 -o smc_analysis/"${mu}"/"$type"/"${pop}_${sub}"/ --spline pchip "$mu" "$smc_in"/"${pop}_${sub}-"*.smc.gz
#smc++ cv --cores 24 -o smc_analysis/"${mu}"/"$type"/"${pop}_${sub}"/ --spline pchip "$mu" "$smc_in"/"${pop}_${sub}-"*.smc.gz

##### FINAL GRAPHICAL OUTPUT #####
echo -e "plot SMC++ results"
date

smc++ plot -c smc_analysis/"$mu"/"$type"/"${pop}_${sub}".pdf smc_analysis/"${mu}"/"$type"/"${pop}_${sub}"/model.final.json
# -c		produces CSV-formatted table containing the data used to generate the plot

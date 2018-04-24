#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemm
#SBATCH -t 24:00:00
#SBATCH -a 26
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=80G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########
## initial set up borrowed from M Stetter ##

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft)
module load conda3
module load vcftools
module load tabix

#############################
### Set up the parameters ###
#############################

####### ARRAY VARIABLES #######
x=$SLURM_ARRAY_TASK_ID
g=$(( x-1 ))


####### PATHS #######
# genome reference file
#genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# Joint VCF file
#vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/all_jointcalls.vcf"
# filtered joint VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# SMC++ prepped file directory
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"

# smc file in $vcf_filt directory
smc_file="smcpp_prunus_biallelic.recode.vcf.gz"

####### PARAMETERS #######
mu="7.77e-9"	# population mutation rate
  # 7.77e-9 (parent to selfed progeny)
  # 9.48e-9 (low heterozygosity peach to progeny)
  # 1.38e-8 (high heterozygosity peach to interspecific cross to selfed progeny)
  # Xie et al. 2016
cut="5000"	# cutoff length for homozygosity

set up i from slurm array ID

####################
### Begin script ###
####################

##### ACQUIRE SETUP DATA #####
# path to sample list
list="/home/dmvelasc/Projects/Prunus/Script/smcpp_data.txt"

# mapfile to extract sample ID and read name information, each line is array item
mapfile -s "$g" -n 1 -t id < "${list}"
# -s number of rows to skip
# -n number of rows to read
# -t (remove leading/trailing whitespace?)
# id is the array name (anything in this position is the array name)

# create an array from each two column line
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop="${arr[0]}"
sub="${arr[1]}"
samples="${arr[2]}"


echo -e "begin file preparation for SMC++\n run for loop by chromosome as per smcpp instructions\n select individuals and populations at this step"
date

##### SMC++ FINAL PREP #####
# *.smc.gz files for each chromosome are the SMC++ output file
for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$smc_file" "$smc_in"/"$pop"_"$sub"-"$i".smc.gz scaffold_"$i" "$pop":"$samples"
done

mkdir -p smc_analysis/"$pop"_"$sub"_"$mu"

echo -e "begin SMC++ analysis"
date

##### SMC++ ANALYSIS #####
# model.final.json and iterations are outputs
smc++ estimate -o smc_analysis/"$pop"_"$sub"_"$mu"/ "$mu" "$smc_in"/"$pop"_"$sub"-*.smc.gz

# --polarization-error 0.5
# --polarization-error: if the identity of the ancestral allele is not known,
# these options can be used to specify a prior over it. With polarization error p,
# emissions probabilities for entry CSFS(a,b) will be computed as
# (1-p) CSFS(a,b) + p CSFS(2-a, n-b). The default setting is 0.5,
# i.e. the identity of the ancestral allele is not known.
# --unfold is an alias for --polarization-error 0. If the ancestral allele is known
# (from an outgroup, say) then this option will use the unfolded SFS for computing
# probabilities. Incorrect usage of this feature may lead to erroneous results.
# $mu is per generation mutation rate, will probably need to run with three different values based on Xie et al.

##### FINAL GRAPHICAL OUTPUT #####
echo -e "plot SMC++ results"
date

smc++ plot -c smc_analysis/"$pop"_"$sub"_"$mu"/"$pop"_"$sub"_"$mu".pdf smc_analysis/"$pop"_"$sub"_"$mu"/model.final.json
smc++ plot -g 10 smc_analysis/"$pop"_"$sub"_"$mu"/"$pop"_"$mu"_"$sub"_years.pdf smc_analysis/"$pop"_"$sub"_"$mu"/model.final.json
#  smc++ plot --logy "$pop"_"$mu"_"$sub"_logY.pdf smc_analysis/model.final.json
# -g		sets generation time in years to scale x-axis, otherwise in coalescent units
# --logy	plots the y-axis on a log scale <- gives an error
# -c		produces CSV-formatted table containing the data used to generate the plot

#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemh
#SBATCH -t 8:00:00
#SBATCH -a 1-40%2
#SBATCH -n 1
#SBATCH -c 14
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=110G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft)
module use /home/dmvelasc/MyModules
module load miniconda3

# Activate Environment
set +u
source /home/dmvelasc/.virtualenvs/smcpp/bin/activate
set -u

########################
### Set up variables ###
########################
x=$SLURM_ARRAY_TASK_ID
g=$(( x-1 ))

####### PATHS #######
# filtered joint VCF file - dulcis test VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# smc file in $vcf_filt directory
smc_file="smcpp_prunus_biallelic.recode.vcf.gz"
# command used [estimate (est) or cross validation (cv)] to create population marginal estimates
type="est"

####### SMC++ PARAMETERS #######
mu="7.77e-9"	# population mutation rate
cut="5000"	# cutoff length for homozygosity
# bed file for masking positions
#mask_file="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.softmasked.bed.bgz"

##############################
##### ACQUIRE SETUP DATA #####
# path to sample list
list="/home/dmvelasc/Projects/Prunus/Script/smcpp_data.txt"
# path to subset pairs for split comparison
pairs="/home/dmvelasc/Projects/Prunus/Script/smcpp_split.txt"


####### ARRAY TO SELECT POPULATIONS FOR COMPARISON #######
echo -e "Selecting populations for split comparison..."
# mapfile to extract first pop/subpop information and send to an array
mapfile -s "$g" -n 1 -t id < "${pairs}"
lines=(`echo "${id[0]}"`)

# original order in smcpp_split.txt
line1="${lines[3]}"	#line number of first pop/subpop wanted
line2="${lines[7]}"	#line number of second pop/subpop wanted

# reverse order in smcpp_split.txt, to test if order matters
#line1="${lines[7]}"	#line number of first pop/subpop wanted
#line2="${lines[3]}"	#line number of second pop/subpop wanted

# original order
echo -e "Selecting populations: ${lines[0]}-${lines[1]} and ${lines[4]}-${lines[5]},\nrepresenting ${lines[2]} and ${lines[6]}, respectively."

# reverse order
#echo -e "Selecting populations: ${lines[4]}-${lines[5]} and ${lines[0]}-${lines[1]},\nrepresenting ${lines[6]} and ${lines[2]}, respectively."
date

####### ARRAY TO EXTRACT POPULATION & SAMPLE INFORMATION FOR SMC++ #######
echo -e "Extract information from population information file"
date

echo -e "Processing population 1"
# mapfile to extract first pop/subpop information and send to an array
mapfile -s $(( line1-1 )) -n 1 -t id < "${list}"
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop1="${arr[0]}"
sub1="${arr[1]}"
samples1="${arr[2]}"

echo -e "Population: $pop1 ; subpopulation: $sub1 ; samples: $samples1"
date

echo -e "Processing population 2"
# mapfile to extract second pop/subpop information and send to an array
mapfile -s $(( line2-1 )) -n 1 -t id < "${list}"
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop2="${arr[0]}"
sub2="${arr[1]}"
samples2="${arr[2]}"

echo -e "Population: $pop1 ; subpopulation: $sub1 ; samples: $samples2"
date
######################

####################
### Begin script ###
####################

mkdir -p smc_analysis/"$mu"/"$type"_cut_split/"$pop1"-"$line1"_"$pop2"-"$line2"
mkdir -p smc_analysis/"$mu"/"$type"_cut_split/"$pop2"-"$line2"_"$pop1"-"$line1"
mkdir -p smc_analysis/"$mu"/"$type"_cut_split/"$pop2"-"$line2"_"$pop1"-"$line1"


################## SPLIT #####################
echo -e "SMC++: create reciprocal joint frequency spectrum datasets and run split analysis"
date

# Create reciprocal joint frequency spectra
for i in {1..8}; do
#  smc++ vcf2smc --mask "$mask_file" "$vcf_filt"/"$smc_file" \
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$smc_file" \
smc_analysis/"$mu"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/"${pop1}_${sub1}-${pop2}_${sub2}-${i}".smc.gz scaffold_"$i" \
"${pop1}:${samples1}" "${pop2}:${samples2}"
#  smc++ vcf2smc --mask "$mask_file" "$vcf_filt"/"$smc_file" \
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$smc_file" \
smc_analysis/"$mu"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/"${pop2}_${sub2}-${pop1}_${sub1}-${i}".smc.gz scaffold_"$i" \
"${pop2}:${samples2}" "${pop1}:${samples1}"
done

echo -e "SMC++: refine the marginal estimates"
date

# Refine the marginal estimates
# order 1
smc++ split -o smc_analysis/"${mu}"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/ \
 smc_analysis/"$mu"/"$type"/cut/"${pop1}_${sub1}"/model.final.json \
 smc_analysis/"$mu"/"$type"/cut/"${pop2}_${sub2}"/model.final.json \
 smc_analysis/"$mu"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/*.smc.gz
# order 2
smc++ split -o smc_analysis/"${mu}"/"$type"_cut_split/"${pop2}-${line2}_${pop1}-${line1}"/ \
 smc_analysis/"$mu"/"$type"/cut/"${pop2}_${sub2}"/model.final.json \
 smc_analysis/"$mu"/"$type"/cut/"${pop1}_${sub1}"/model.final.json \
 smc_analysis/"$mu"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/*.smc.gz

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
# Order 1
smc++ plot -c \
smc_analysis/"$mu"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/"${pop1}-${line1}_${pop2}-${line2}_${mu}".pdf \
smc_analysis/"$mu"/"$type"_cut_split/"${pop1}-${line1}_${pop2}-${line2}"/model.final.json
# Order 2
smc++ plot -c \
smc_analysis/"$mu"/"$type"_cut_split/"${pop2}-${line2}_${pop1}-${line1}"/"${pop2}-${line2}_${pop1}-${line1}_${mu}".pdf \
smc_analysis/"$mu"/"$type"_cut_split/"${pop2}-${line2}_${pop1}-${line1}"/model.final.json

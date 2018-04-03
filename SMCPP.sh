#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemm
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -c 10
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

####### PATHS #######
# genome reference file
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# Joint VCF file - test joint VCF file
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/all_jointcalls.vcf"
# filtered joint VCF file - dulcis test VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# SMC++ prepped file directory
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"

####### PARAMETERS #######
mu="7.77e-9"	# population mutation rate
  # 7.77e-9 (parent to selfed progeny)
  # 9.48e-9 (low heterozygosity peach to progeny)
  # 1.38e-8 (high heterozygosity peach to interspecific cross to selfed progeny)
  # Xie et al. 2016
cut="5000"	# cutoff length for homozygosity
pop="PV"	# population

####################
### Begin script ###
####################
for i in {0..21}; do

  ##### ACQUIRE SETUP DATA #####
  # path to sample list
  list="/home/dmvelasc/Projects/Prunus/Script/smcpp_data.txt"
  # mapfile to extract sample ID and read name information, each line is array item
  mapfile -s "$i" -n 1 -t id < "${list}"
  # -s number of rows to skip
  # -n number of rows to read
  # -t (remove leading/trailing whitespace?)
  # id is the array name (anything in this position is the array name if nothing then called ?array)
  # create an array from each two column line
  arr=(`echo "${id[0]}"`)

  # declare variables, created from array
  pop="${arr[0]}"
  sub="${arr[1]}"
  sample1="${arr[2]}"
  sample2="${arr[3]}"
  sample3="${arr[4]}"
  sample4="${arr[5]}"


  echo -e "begin SMC++ preparation\n get individuals"
  date

  ##### NEEDED FOR INITIAL PREP #####

  # PD; dulcis; subset=all; 18 individuals; 12 CPU
  # PP; persica; subset=all; 14 individuals; 12 CPU
  # PM; mira; subset=all; 6 individuals; 6 CPU
  # PV; davidiana; subset=all; 6 individuals; 6 CPU; 4 individuals 10 CPU (????)
  # PS; kansuensis; subset=all; 4 individuals; 5 CPU
  # PG; ferganensis; subset=all; 4 individuals; 4 CPU
  vcftools --vcf "$vcf" --indv "$sample1" --indv "$sample2" --indv "$sample3" --indv "$sample4" --min-alleles 2 --max-alleles 2 --recode --out "$sub"_"$pop"

  ##### NEEDED FOR INITIAL PREP #####
  mv /home/dmvelasc/Projects/Prunus/Analysis/smcpp/"$sub"_"$pop".recode.vcf "$vcf_filt"/

  echo -e "convert vcf file to SMC++ format file"
  date

  ##### NEEDED FOR INITIAL PREP #####
  bgzip -f "$vcf_filt"/"$sub"_"$pop".recode.vcf > "$vcf_filt"/"$sub"_"$pop".recode.vcf.gz
  tabix -fp vcf "$vcf_filt"/"$sub"_"$pop".recode.vcf.gz

  echo -e "Run for loop by chromosome as per smcpp instructions"
  date

  ##### NEEDE FOR INITIAL PREP #####
  for i in {1..8}; do
    smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/"$sub"_"$pop".recode.vcf.gz "$smc_in"/"$sub"_"$pop"_"$i".smc.gz scaffold_"$i" "$pop":"$sample1","$sample2","$sample3","$sample4"
  done


  echo -e "begin SMC++ analysis"
  date
  ##### NEEDED FOR ANALYSIS #####
  # SMC++ analysis
  smc++ estimate -o smc_analysis/ "$mu" "$smc_in"/*"$pop"*.smc.gz

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
  smc++ plot -c "$pop"_"$mu"_"$sub".pdf smc_analysis/model.final.json
  smc++ plot -g 10 "$pop"_"$mu"_"$sub"_years.pdf smc_analysis/model.final.json
  smc++ plot --logy "$pop"_"$mu"_"$sub"_logY.pdf smc_analysis/model.final.json
  #-g	sets generation time in years to scale x-axis, otherwise in coalescent units
  #--logy	plots the y-axis on a log scale
  #-c	produces CSV-formatted table containing the data used to generate the plot

  # move files output files to subdirectory
  mkdir -p smc_analysis/"$pop"_"$sub"_"$mu"
  mv smc_analysis/model.final.json smc_analysis/"$pop"_"$sub"_"$mu"
  mv smc_analysis/.model.iter*.json smc_analysis/"$pop"_"$sub"_"$mu"
  mv smc_analysis/.debug.txt smc_analysis/"$pop"_"$sub"_"$mu"

done

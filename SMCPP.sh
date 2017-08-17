#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemh
#SBATCH -t 4:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=56G
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
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/all_dulcis.recode.vcf"
# SMC++ prepped file
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"
final="mu1_38x-8_alldulcis" #final directory name

####### PARAMETERS #######
mu="1.38e-8"	# population mutation rate
cut="5000"	# cutoff length for homozygosity
pop="dulcis"	# population

####################
### Begin script ###
####################
echo -e "begin SMC++ preparation\n get individuals"
date

# select individuals
vcftools --vcf "$vcf" --indv PD02 --indv PD03 --indv PD04 --indv PD05 --indv PD06 --indv PD07 --indv PD08 --indv PD09 --indv PD10 --indv PD11 --indv PD12 --indv PD13 --indv PD14 --indv PD16 --indv PD17 --indv PD18 --indv PD20 --indv PD21 --min-alleles 2 --max-alleles 2 --recode --out all_"$pop"
mv /home/dmvelasc/Projects/Prunus/Analysis/smcpp/all_"$pop".recode.vcf /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/

echo -e "convert vcf file to SMC++ format file"
date

bgzip -f "$vcf_filt" > "$vcf_filt".gz
tabix -fp vcf "$vcf_filt".gz

echo -e "Run for loop by chromosome as per smcpp instructions"
date

for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt".gz "$smc_in"/all_"$pop"_"$i".smc.gz scaffold_"$i" "$pop":PD02,PD03,PD04,PD05,PD06,PD07,PD08,PD09,PD10,PD11,PD12,PD13,PD14,PD16,PD17,PD18,PD20,PD21
done


echo -e "begin SMC++ analysis"
date
# SMC++ analysis
smc++ estimate -o smc_analysis/ "$mu" "$smc_in"/*.smc.gz

#--polarization-error 0.5
# --polarization-error: if the identity of the ancestral allele is not known,
# these options can be used to specify a prior over it. With polarization error p,
# emissions probabilities for entry CSFS(a,b) will be computed as
# (1-p) CSFS(a,b) + p CSFS(2-a, n-b). The default setting is 0.5,
# i.e. the identity of the ancestral allele is not known.
# --unfold is an alias for --polarization-error 0. If the ancestral allele is known
# (from an outgroup, say) then this option will use the unfolded SFS for computing
# probabilities. Incorrect usage of this feature may lead to erroneous results.
# $mu is per generation mutation rate, will probably need to run with three different values based on Xie et al.

echo -e "plot SMC++ results"
date
smc++ plot history_"$pop"_"$mu".pdf smc_analysis/model.final.json
#-g	sets generation time in years to scale x-axis, otherwise in coalescent units
#--logy	plots the y-axis on a log scale
#-c	produces CSV-formatted table containing the data used to generate the plot

# move files output files to subdirectory
mkdir -p smc_analysis/"$final"
mv smc_analysis/model.final.json smc_analysis/"$final"
mv smc_analysis/.model.iter*.json smc_analysis/"$final"
mv smc_analysis/.debug.txt smc_analysis/"$final"

#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemm
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=92G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

# Prepare the VCF file(s) for SMC++ analysis

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

####### PARAMETERS #######

####################
### Begin script ###
####################
echo -e "begin SMC++ preparation\n get individuals"
date

# select individuals
vcftools --vcf "$vcf" --indv PD02 --indv PD03 --indv PD04 --indv PD05 --indv PD06 --indv PD07 --indv PD08 --indv PD09 --indv PD10 --indv PD11 --indv PD12 --indv PD13 --indv PD14 --indv PD16 --indv PD17 --indv PD18 --indv PD20 --indv PD21 --indv PP02 --indv PP03 --indv PP04 --indv PP05 --indv PP06 --indv PP08 --indv PP11 --indv PP13 --indv PP14 --indv PP15 --indv PP37 --indv PP38 --indv PP39 --indv PP40 --indv PM01 --indv PM02 --indv PM03 --indv PM04 --indv PM05 --indv PM06 --indv PV01 --indv PV02 --indv PV03 --indv PV04 --indv PV05 --indv PV06 --indv PS01 --indv PS02 --indv PS03 --indv PS04 --indv PG02 --indv PG03 --indv PG04 --indv PG05 --min-alleles 2 --max-alleles 2 --recode --out smcpp_prunus_biallelic

mv /home/dmvelasc/Projects/Prunus/Analysis/smcpp/smcpp_prunus_biallelic.recode.vcf "$vcf_filt"/

echo -e "convert vcf file to SMC++ format file"
date

bgzip -f "$vcf_filt"/smcpp_prunus_biallelic.recode.vcf > "$vcf_filt"/smcpp_prunus_biallelic.recode.vcf.gz
tabix -fp vcf "$vcf_filt"/smcpp_prunus_biallelic.recode.vcf.gz

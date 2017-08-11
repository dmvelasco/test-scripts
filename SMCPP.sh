#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemm
#SBATCH -t 10-00:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=56G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########
###### borrowed heavily from M Stetter #####

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
# Declare number variables
#x=$SLURM_ARRAY_TASK_ID
#i=$(( x-1 ))

# full set of initial samples
#declare -a id=(PC01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PR01 PS02)
#sample="${id["$i"]}"

# path for genome reference file
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# path for joint VCF file - test joint VCF file
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_jointcalls.vcf"
# path for filtered joint VCF file - dulcis test VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_dulcis.recode.vcf"
# path for SMC++ prepped file
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"


# Declare directories and file prefix
#dir1="/home/dmvelasc/Software/bedtools2/bin"		# bedtools location
#dir2="/home/dmvelasc/Projects/Prunus/Data/BAM"		# BAM directory
dir3="/home/dmvelasc/Data/references/persica-SCF"	# gff3 location
dir4="/scratch/dmvelasc"

####################
### Begin script ###
####################
echo -e "begin SMC++ preparation\n get individuals"
date

# select individuals
vcftools --vcf "$vcf" --indv PD03 --indv PD04 --indv PD05 --indv PD06 --indv PD07 --indv PD08 --indv PD09 --min-alleles 2 --max-alleles 2 --recode --out test_dulcis
mv /home/dmvelasc/Projects/Prunus/Analysis/smcpp/test_dulcis.recode.vcf /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/

echo -e "convert vcf file to SMC++ format file"
date

bgzip -f "$vcf_filt" > "$vcf_filt".gz
tabix -fp vcf "$vcf_filt".gz
#smc++ vcf2smc --missing-cutoff 1000 --length 80000000 "$vcf_filt".gz "$smc_in"/test_dulcis_1.smc.gz scaffold_1 Pop1:PD03,PD04,PD05,PD06,PD07,PD08,PD09
# above is modified vcf2smc line, however, put into for loop below to iterate over first eight scaffolds of Prunus

# Markus' original SMC++ line
#smc++ vcf2smc --missing-cutoff 1000 --length 80000000 chr10_282_4ind.vcf.gz smc_input/chr10_282_4ind.smc.gz 10 Pop1:282set_B73,282set_W22,282set_Ki11,282set_T8
#after smc_input/file 10 Pop1:232set_B73,282set_W22,282set_Ki11,282set_T8
# 10	chromosome number
# Pop1:282set_B73,282set_W22,282set_Ki11,282set_T8
# population and members of the population
# --missing-cutoff is the length of runs of homozygosity to ignore, assumes missing data
# --length is length of the chromosome, default uses header length <- does vcf tools or other manipulations remove header?

echo -e "Run for loop by chromosome as per smcpp instructions"
date

for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff 5000 "$vcf_filt".gz "$smc_in"/test_dulcis_"$i".smc.gz scaffold_"$i" Pop1:PD03,PD04,PD05,PD06,PD07,PD08,PD09
done


echo -e "begin SMC++ analysis"
date
# SMC++ analysis
smc++ estimate -o smc_analysis/ 3e-8 "$smc_in"/*.smc.gz
#--polarization-error 0.5
# --polarization-error: if the identity of the ancestral allele is not known,
# these options can be used to specify a prior over it. With polarization error p,
# emissions probabilities for entry CSFS(a,b) will be computed as
# (1-p) CSFS(a,b) + p CSFS(2-a, n-b). The default setting is 0.5,
# i.e. the identity of the ancestral allele is not known.
# --unfold is an alias for --polarization-error 0. If the ancestral allele is known
# (from an outgroup, say) then this option will use the unfolded SFS for computing
# probabilities. Incorrect usage of this feature may lead to erroneous results.
#3e-8 is per generation mutation rate, will probably need to run with three different values based on Xie et al.

echo -e "plot SMC++ results"
date
smc++ plot history_dulcis.pdf smc_analysis/model.final.json
#

#-g	sets generation time in years to scale x-axis, otherwise in coalescent units
#--logy	plots the y-axis on a log scale
#-c	produces CSV-formatted table containing the data used to generate the plot



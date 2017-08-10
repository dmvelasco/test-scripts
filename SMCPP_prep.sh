#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemm
#SBATCH -t 10-00:00:00
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=8G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft)
module load conda3
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

# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# file directory for $sample vcf to convert to fasta
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_jointcalls.vcf"
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/test_jointcalls.smc.gz"


# Declare directories and file prefix
dir1="/home/dmvelasc/Software/bedtools2/bin"		# bedtools location
dir2="/home/dmvelasc/Projects/Prunus/Data/BAM"		# BAM directory
dir3="/home/dmvelasc/Data/references/persica-SCF"	# gff3 location
dir4="/scratch/dmvelasc"

####################
### Begin script ###
####################
echo -e "begin SMC++ preparation\n get individuals"
date

# select individuals
module load vcftools
vcftools --gzvcf "$vcf".gz --indv PD03 --indv PD04 --indv PD05 --indv PD06 --indv PD07 --indv PD08 --indv PD09 --min-alleles 2 --max-alleles 2 --recode --out dulcis

echo -e "begin SMC++ preparation\n convert vcf file to SMC++ format file"
date

bgzip "$vcf" > "$vcf".gz
tabix -fp vcf "$vcf".gz
smc++ vcf2smc --missing-cutoff 1000 --length 80000000 "$vcf".gz "$smc_in" scaffold_1 Pop1:PD03,PD04,PD05,PD06,PD07,PD08,PD09
# Markus' original SMC++ line
#smc++ vcf2smc --missing-cutoff 1000 --length 80000000 chr10_282_4ind.vcf.gz smc_input/chr10_282_4ind.smc.gz 10 Pop1:282set_B73,282set_W22,282set_Ki11,282set_T8
#after smc file 10 Pop1:232set_B73,282set_W22,282set_Ki11,282set_T8
# 10	chromosome number
# Pop1:282set_B73,282set_W22,282set_Ki11,282set_T8
# population and members of the population


#run array with each chromosome, not sure needed for prunus b/c small genome
#SBATCH -a 1-10

echo -e "begin SMC++ preparation\n for each chromosome"
date

scaf=$SLURM_ARRAY_TASK_ID
module load tabix

bgzip -c > "$vcf".gz
tabix -fp vcf "$vcf".gz


echo -e "begin SMC++ analysis"
date
# SMC++ analysis
smc++ estimate -o smc_analysis/ 3e-8 smc_input/*.smc.gz
#--polarization-error 0.5
# --polarization-error: if the identity of the ancestral allele is not known,
# these options can be used to specify a prior over it. With polarization error p,
# emissions probabilities for entry CSFS(a,b) will be computed as
# (1-p) CSFS(a,b) + p CSFS(2-a, n-b). The default setting is 0.5,
# i.e. the identity of the ancestral allele is not known.
# --unfold is an alias for --polarization-error 0. If the ancestral allele is known
# (from an outgroup, say) then this option will use the unfolded SFS for computing
# probabilities. Incorrect usage of this feature may lead to erroneous results.
#3e-8 is per generation mutation rate, will probably need to run with three different values (Xie et al.)

echo -e "plot SMC++ results"
date
smc++ plot history_all_282.pdf analysis/model.final.json
#

#-g	sets generation time in years to scale x-axis, otherwise in coalescent units
#--logy	plots the y-axis on a log scale
#-c	produces CSV-formatted table containing the data used to generate the plot

###### From M Stetter

###### Convert VCF
#!/bin/bash -l
#SBATCH -D /home/mstetter/smcpp_test/
#SBATCH -o /home/mstetter/smcpp_test/stdout-%j.txt
#SBATCH -J extract_282
#SBATCH -t 00:05:00
set -e
set -u

module load tabix

source /home/mstetter/tools/smcpp/bin/activate
bgzip chr10_282_4ind.vcf
tabix -fp vcf chr10_282_4ind.vcf.gz
smc++ vcf2smc --missing-cutoff 1000 --length 80000000 chr10_282_4ind.vcf.gz smc_input/chr10_282_4ind.smc.gz 10 Pop1:282set_B73,282set_W22,282set_Ki11,282set_T8



#!/bin/bash -l
#SBATCH -D /home/mstetter/282_analysis/
#SBATCH -o /home/mstetter/logs/smcpp-%j.txt
#SBATCH -J vcf2smc
#SBATCH -t 2-00:00
#SBATCH --array=1-10
set -e
set -u


hset=$SLURM_ARRAY_TASK_ID
module load tabix

gunzip -c /group/jrigrp/Share/genotypes/282_7X/c${hset}_282_corrected_onHmp321.vcf.gz | bgzip -c > c${hset}_282_corrected_onHmp321.vcf.gz
tabix -fp vcf c${hset}_282_corrected_onHmp321.vcf.gz


###### Get individuals
#!/bin/bash -l
#SBATCH -D /home/mstetter/smcpp_test/
#SBATCH -o /home/mstetter/smcpp_test/stdout-%j.txt
#SBATCH -J extract_282
#SBATCH -t 00:05:00
set -e
set -u

module load vcftools
vcftools --gzvcf /group/jrigrp/Share/genotypes/282_7X/c10_282_corrected_onHmp321.vcf.gz --indv 282set_B73 --indv 282set_W22 --indv 282set_Ki11 --indv 282set_T8 --min-alleles 2 --max-alleles 2 --recode --out chr10_282_4ind

##### SMC analysis
#!/bin/bash -l
#SBATCH -D /home/mstetter/smcpp_test/
#SBATCH -o /home/mstetter/smcpp_test/stdout-%j.txt
#SBATCH -J smcpp_analysis
#SBATCH -t 04:00:00
#SBATCH -p bigmemm
set -e
set -u

source /home/mstetter/tools/smcpp/bin/activate

smc++ estimate -o smc_analysis/ --polarization-error 0.5 3e-8 smc_input/*.smc.gz
smc++ plot history_all_282.pdf analysis/model.final.json

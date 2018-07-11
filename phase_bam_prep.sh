#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcffilt.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcffilt.txt
#SBATCH -J fasta
#SBATCH -p bigmemh
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=16G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

####################
### Load modules ###
####################
module load vcftools
module load tabix

####################
### Begin script ###
####################
echo -e "recode jointcalls to limit to 2 alleles, should eliminate asterisk alleles, which indicate a spanning deletion"
date

# select individuals
vcftools --vcf all_jointcalls.vcf --max-alleles 2 --recode --out all_jointcalls_biallelic
bgzip -f all_jointcalls_biallelic.recode.vcf > all_jointcalls_biallelic.recode.vcf.gz
tabix -fp vcf all_jointcalls_biallelic.recode.vcf.gz

#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-phase.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-phase.txt
#SBATCH -J test
#SBATCH -p bigmemm
#SBATCH -t 14-00:00:00
#SBATCH -n 1
#SBATCH -c 38
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=296G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft and WhatsHap)
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
bam_dir="/group/jrigrp3/Velasco/Prunus/BAM"

####################
### Begin script ###
####################

# basic script
srun whatshap phase -o "$vcf_filt"/single_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
"$bam_dir"/PB01_sorted_markdup.bam \
"$bam_dir"/PR01_sorted_markdup.bam \
"$bam_dir"/PU01_sorted_markdup.bam \
"$bam_dir"/PC01_sorted_markdup.bam \
"$bam_dir"/PF01_sorted_markdup.bam \
"$bam_dir"/PK01_sorted_markdup.bam \
"$bam_dir"/PT01_sorted_markdup.bam \
"$bam_dir"/PV01_sorted_markdup.bam \
"$bam_dir"/PV02_sorted_markdup.bam \
"$bam_dir"/PV03_sorted_markdup.bam \
"$bam_dir"/PV04_sorted_markdup.bam \
"$bam_dir"/PV05_sorted_markdup.bam \
"$bam_dir"/PV06_sorted_markdup.bam \
"$bam_dir"/PG02_sorted_markdup.bam \
"$bam_dir"/PG03_sorted_markdup.bam \
"$bam_dir"/PG04_sorted_markdup.bam \
"$bam_dir"/PG05_sorted_markdup.bam \
"$bam_dir"/PS01_sorted_markdup.bam \
"$bam_dir"/PS02_sorted_markdup.bam \
"$bam_dir"/PS03_sorted_markdup.bam \
"$bam_dir"/PS04_sorted_markdup.bam \
"$bam_dir"/PM01_sorted_markdup.bam \
"$bam_dir"/PM02_sorted_markdup.bam \
"$bam_dir"/PM03_sorted_markdup.bam \
"$bam_dir"/PM04_sorted_markdup.bam \
"$bam_dir"/PM05_sorted_markdup.bam \
"$bam_dir"/PM06_sorted_markdup.bam \
"$bam_dir"/PP01_sorted_markdup.bam \
"$bam_dir"/PP02_sorted_markdup.bam \
"$bam_dir"/PP03_sorted_markdup.bam \
"$bam_dir"/PP04_sorted_markdup.bam \
"$bam_dir"/PP05_sorted_markdup.bam \
"$bam_dir"/PP06_sorted_markdup.bam \
"$bam_dir"/PP08_sorted_markdup.bam \
"$bam_dir"/PP09_sorted_markdup.bam \
"$bam_dir"/PP11_sorted_markdup.bam \
"$bam_dir"/PP12_sorted_markdup.bam \
"$bam_dir"/PP13_sorted_markdup.bam \
"$bam_dir"/PP14_sorted_markdup.bam \
"$bam_dir"/PP15_sorted_markdup.bam \
"$bam_dir"/PP37_sorted_markdup.bam \
"$bam_dir"/PP38_sorted_markdup.bam \
"$bam_dir"/PP39_sorted_markdup.bam \
"$bam_dir"/PP40_sorted_markdup.bam \
"$bam_dir"/PD01_sorted_markdup.bam \
"$bam_dir"/PD02_sorted_markdup.bam \
"$bam_dir"/PD03_sorted_markdup.bam \
"$bam_dir"/PD04_sorted_markdup.bam \
"$bam_dir"/PD05_sorted_markdup.bam \
"$bam_dir"/PD06_sorted_markdup.bam \
"$bam_dir"/PD07_sorted_markdup.bam \
"$bam_dir"/PD08_sorted_markdup.bam \
"$bam_dir"/PD09_sorted_markdup.bam \
"$bam_dir"/PD10_sorted_markdup.bam \
"$bam_dir"/PD11_sorted_markdup.bam \
"$bam_dir"/PD12_sorted_markdup.bam \
"$bam_dir"/PD13_sorted_markdup.bam \
"$bam_dir"/PD14_sorted_markdup.bam \
"$bam_dir"/PD16_sorted_markdup.bam \
"$bam_dir"/PD17_sorted_markdup.bam \
"$bam_dir"/PD18_sorted_markdup.bam \
"$bam_dir"/PD19_sorted_markdup.bam \
"$bam_dir"/PD20_sorted_markdup.bam \
"$bam_dir"/PD21_sorted_markdup.bam \
--indels

# multiple bams can be combined on the command line, program automatically detects sample(s) from file
#Input pre-procession, selection, and filtering:
#  --max-coverage MAXCOV, -H MAXCOV
#                        Reduce coverage to at most MAXCOV (default: 15).
#  --mapping-quality QUAL, --mapq QUAL
#                        Minimum mapping quality (default: 20)
#  --indels              Also phase indels (default: do not phase indels)
#  --ignore-read-groups  Ignore read groups in BAM header and assume all reads
#                        come from the same sample.
#  --sample SAMPLE       Name of a sample to phase. If not given, all samples
#                        in the input VCF are phased. Can be used multiple
#                        times.
#Genotyping:
#  --full-genotyping     Completely re-genotype all variants based on read
#                        data, ignores all genotype data that might be present
#                        in the VCF (EXPERIMENTAL FEATURE).
#  --distrust-genotypes  Allow switching variants from hetero- to homozygous in
#                        an optimal solution (see documentation).
#  --include-homozygous  Also work on homozygous variants, which might be
#                        turned to heterozygous
#  --changed-genotype-list FILE
#                        Write list of changed genotypes to FILE.


#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-phase.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-phase.txt
#SBATCH -J test
#SBATCH -p bigmemm
#SBATCH -t 14-00:00:00
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=152G
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
--sample PB01 \
--sample PR01 \
--sample PU01 \
--sample PC01 \
--sample PF01 \
--sample PK01 \
--sample PT01 \
--indels

srun whatshap phase -o "$vcf_filt"/PV_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
"$bam_dir"/PV01_sorted_markdup.bam \
"$bam_dir"/PV02_sorted_markdup.bam \
"$bam_dir"/PV03_sorted_markdup.bam \
"$bam_dir"/PV04_sorted_markdup.bam \
"$bam_dir"/PV05_sorted_markdup.bam \
"$bam_dir"/PV06_sorted_markdup.bam \
--sample PV01 \
--sample PV02 \
--sample PV03 \
--sample PV04 \
--sample PV05 \
--sample PV06 \
--indels

srun whatshap phase -o "$vcf_filt"/PG_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
"$bam_dir"/PG02_sorted_markdup.bam \
"$bam_dir"/PG03_sorted_markdup.bam \
"$bam_dir"/PG04_sorted_markdup.bam \
"$bam_dir"/PG05_sorted_markdup.bam \
--sample PG02 \
--sample PG03 \
--sample PG04 \
--sample PG05 \
--indels

srun whatshap phase -o "$vcf_filt"/PS_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
--sample PS01 \
--sample PS02 \
--sample PS03 \
--sample PS04 \
"$bam_dir"/PS01_sorted_markdup.bam \
"$bam_dir"/PS02_sorted_markdup.bam \
"$bam_dir"/PS03_sorted_markdup.bam \
"$bam_dir"/PS04_sorted_markdup.bam \
--indels

srun whatshap phase -o "$vcf_filt"/PM_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
"$bam_dir"/PM01_sorted_markdup.bam \
"$bam_dir"/PM02_sorted_markdup.bam \
"$bam_dir"/PM03_sorted_markdup.bam \
"$bam_dir"/PM04_sorted_markdup.bam \
"$bam_dir"/PM05_sorted_markdup.bam \
"$bam_dir"/PM06_sorted_markdup.bam \
--sample PM01 \
--sample PM02 \
--sample PM03 \
--sample PM04 \
--sample PM05 \
--sample PM06 \
--indels

srun whatshap phase -o "$vcf_filt"/PP_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
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
--sample PP01 \
--sample PP02 \
--sample PP03 \
--sample PP04 \
--sample PP05 \
--sample PP06 \
--sample PP08 \
--sample PP09 \
--sample PP11 \
--sample PP12 \
--sample PP13 \
--sample PP14 \
--sample PP15 \
--sample PP37 \
--sample PP38 \
--sample PP39 \
--sample PP40 \
--indels

srun whatshap phase -o "$vcf_filt"/PD_jointcalls_phased.vcf.gz \
-r "$genome" "$vcf" \
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
--sample PD01 \
--sample PD02 \
--sample PD03 \
--sample PD04 \
--sample PD05 \
--sample PD06 \
--sample PD07 \
--sample PD08 \
--sample PD09 \
--sample PD10 \
--sample PD11 \
--sample PD12 \
--sample PD13 \
--sample PD14 \
--sample PD16 \
--sample PD17 \
--sample PD18 \
--sample PD19 \
--sample PD20 \
--sample PD21 \
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


#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis/bcf
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf2ped.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf2ped.txt
#SBATCH -J convert
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u


# %j is job allocation number

# Declare directories
bin="/home/dmvelasc/bin"

# import genotype data from vcf to plink
"$bin"/plink --vcf gbs_plink.vcf --const-fid 0 --biallelic-only strict list --vcf-min-qual 30 --allow-extra-chr --chr-set 34 --allow-no-sex --make-bed --out convert

# --const-fid converts sample IDs to individula IDs while setting all family IDs to single value (default 0)
# --biallelic-only keeps only biallelic loci
#                  strict - removes any loci with more than one alternate even if no genotypes present with 3rd or greater allele
#                  list - prints a list of loci not kept
# --vcf-min-qual utilizes the score in the QUAL column not the FILTER column, which is used for --vcf-filter
#                not needed because gbs vcf file all qual are greater than 20
# --recode directs it to create .ped file set (.ped, .map)
#            can interchange with alternate declarations
# --out indicates outfile prefix
#       convert is file prefix for new files
# --allow-extra-chromosome allows more chromosomes than in chromosome set
# --chr-set chromosomes higher than this number are disallowed and program exits
# --allow-no-sex allows ambiguous or non-identified sex


# remove loci with >10% missing data and genotypes with >75% missing data
"$bin"/plink --bfile convert --chr-set 34 --geno 0.1 --allow-no-sex --make-bed --out qc

# LD pruning
"$bin"/plink --bfile qc --allow-no-sex --allow-extra-chr --chr-set 34 --indep-pairwise 10 50 0.4

# original parameters 1000 50 0.05 => window in Kb, window, r2 <= window of 1 Mb too big, 10 Kb more likely


# creates bed files from pruned data
"$bin"/plink --bfile qc --chr-set 34 --extract plink.prune.in --allow-no-sex --recode --out ld_pruned

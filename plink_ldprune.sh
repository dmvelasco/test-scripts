#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-plink.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-plink.txt
#SBATCH -J assoc
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2:00:00
# -t time hours:min:sec
# how are days set? maximum 180 days = 180days*24hrs/day=4320hrs
set -e
set -u


# Load zlib
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"				# program directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"      # VCF directory
infile="prunus-A.flt_test.vcf.bzip"			# original file
out1="prunus_plink-2"					# first outfile
out2="plink_test-2"					# second outfile
pruned="ld_pruned-10"					# final pruned file

# filter vcf for snps only
"$dir1"/bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 -i "%QUAL>=30 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out1".vcf -v snps "$dir2"/"$infile"
# -i, --include EXPRESSION
# above for variant quality, depth, and genotype quality
# -r, --regions chr|chr:pos|chr:from-to|chr:from-[,…]
# Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.

# remove "scaffold_"
perl -plne 's/scaffold_(\w+)/$1/' "$out1".vcf > "$out2".vcf


# import genotype data from vcf to plink
"$dir1"/plink --vcf "$out2".vcf --const-fid 0 --chr-set 8 --allow-no-sex --make-bed --out convert

# --const-fid converts sample IDs to individula IDs while setting all family IDs to single value (default 0)
# --make-bed directs it to create .bed file set (.bed, .bim, .fam)
#            can interchange with alternate declarations
# --out indicates outfile prefix
#       convert is file prefix for new files
# --chr-set chromosomes higher than this number are disallowed and program exits
# --allow-no-sex allows ambiguous or non-identified sex


# LD pruning
# --indep-pairwise [window size]<kb> [step size (locus ct)] [r^2 threshold]
# window size in variant count or kilobase (if the 'kb' modifier is present)  OOPS, was not using kb flag
#first test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 10 5 0.4
#second test ld_pruned-1
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 3 5 0.4
#third test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 10 10 0.4
#fourth test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 10 5 0.2
#fifth test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 1 5 0.4
#sixth test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 10 10 0.2
#seventh test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 2 5 0.4
#eigth test ld_pruned
#"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 2 10 0.4
#ninth test ld_pruned
"$dir1"/plink --bfile convert --allow-no-sex --chr-set 8 --indep-pairwise 2kb 5 0


# creates bed files from pruned data
#"$dir1"/plink --bfile convert --chr-set 8 --extract plink.prune.in --allow-no-sex --make-bed --out ld_pruned
# creates vcf files from pruned data
"$dir1"/plink --bfile convert --chr-set 8 --extract plink.prune.in --allow-no-sex --recode vcf --out "$pruned"

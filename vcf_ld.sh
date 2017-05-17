#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-thin.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-thin.txt
#SBATCH -J thin
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
out1="first"						# first outfile
window="nowindow"					# size of window for ld calculation, or nowindow for none
out2="noscaffold"					# second outfile
thin="500"						# minimum snp spacing
out3="thin-500"						# final file prefix

# filter vcf for snps only
"$dir1"/bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 -i "%QUAL>=250 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out1".vcf -v snps "$dir2"/"$infile"
# -i, --include EXPRESSION
# above for variant quality, depth, and genotype quality
# -r, --regions chr|chr:pos|chr:from-to|chr:from-[,…]
# Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.

# remove "scaffold_"
perl -plne 's/scaffold_(\w+)/$1/' "$out1".vcf > "$out2".vcf

# thin to at least 250 bp
"$dir1"/vcftools --vcf "$out2".vcf --thin "$thin" --recode --out "$out3"

# get ld statistics with vcf tools, possibly filter utilizing data
#"$dir1"/vcftools --vcf "$out3".vcf --geno-r2 --ld-window-bp "$window" --out "$out3"_"$window"
#"$dir1"/vcftools --vcf "$out3".recode.vcf --geno-r2 --out "$out3"_"$window"

#grep, if needed
#grep -n "0\.02" ld_window_2000.geno.ld > ld_window_2000_low.geno.ld


# LD pruning


# create pruned VCF
#"$pruned"

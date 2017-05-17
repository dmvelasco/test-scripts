#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf2struct.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf2struct.txt
#SBATCH -J vcf2struct
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Load zlib 1.2.8 not recognized
#module load zlib

# Declare directories and file prefix
dir1="/home/dmvelasc/bin"				# software binary directory
dir2="/home/dmvelasc/Software"				# general software directory
file="prunus"						# file prefix

#XXXXXXXXXXXXXXXXXXXXXXXXXXXX

##### FILTER VARIANTS #####
# filter raw VCF SNPs for overall quality of 30, Average depth > 3 and <= 100, GQ >= 20 (approximate average depth is 30X)
"$dir1"/bcftools view -O z -o "$file".GQ20flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %AVG(DP)<=100 && %MIN(GQ)>=20 && INFO/MQ>=30" -q 0.05:minor -v snps "$file".raw.vcf.bzip
# -c/C, --min-ac/--max-ac <int>[:<type>]      minimum/maximum count for non-reference (nref), 1st alternate (alt1) or minor (minor) alleles [nref]
# -q/Q, --min-af/--max-af <float>[:<type>]    minimum/maximum frequency for non-reference (nref), 1st alternate (alt1) or minor (minor) alleles [nref]
# from http://wiki.bits.vib.be/index.php/NGS_Exercise.5#call_variants_with_samtools_and_samtools_bcftools
# include additional filter for MAF would add to -i statement MAF[0]>0.05
# MAF[0]<0.05 means select rare variants at 5% cutoff, think using MAF[0]>0.05 means select variants above 5% cutoff

# changed from GQ>=30 to GQ >=20

#XXXXXXXXXXXXXXXXXXXXXXXXXXXX

##### Add SNP ID to filtered VCF file #####
# parse a VCF file without SNP ID and add; without ID vcf2hmp.pl script does not work properly
# Prunus VCF files using peach v1.0 genome, including scaffolds, remove word scaffold and add SNP ID
# VCF columns of interest
# position:     col1 = $1               col2 = $2               col3 = $3
# contents:     scaffold_#      position        ID (all .)
# desired:      #               position        s#_position

# perl
# 	remove word scaffold in file

# awk
# 	designate file separator as tab
# 	if the first field begins with a '#' print the line
# 	otherwise if the third field is a '.' then make it empty
# 	create a variable holding the concatenation of 's', the first field (chromosome/scaffold number),
# 	'_' and the second field (position) then replace the original '.' in field 3 with
# 	the concatenated value and print the line
# 	OFS ="\t" <- to designate that output file is tabl separated

zcat "$file".GQ20flt.vcf.bzip | perl -plne 's/scaffold_(\w+)/$1/' - | awk '{ OFS="\t"; if ($1 ~ /^[#]/) {print} else if ($3 ~ /\./) { sub(/\$3/, ""); temp = ("s"$1"_"$2); $3 = temp; print } ; }' - > "$file".GQ20fltID.vcf

# source of the complex action for the else if statement is:
# http://asdirkseesit.blogspot.com/2012/09/awk-substitute-and-concatenate-columns.html


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

##### FILTER VCF, CONVERT TO HAPMAP, RANDOMLY SELECT 10000 LOCI, CONVERT TO STRUCTURE FORMAT #####

# Filter VCF and convert to hapmap format file
cat "$file".GQ20fltID.vcf | "$dir2"/vcf2hp2.pl > "$file".GQ20.hmp

# Randomly select 10K SNP loci
grep -v "^#" "$file".GQ20.hmp | sort -R - | tail -10000 > "$file"10K_nohdr.GQ20.hmp

# Restore header to file
grep "^rs#" "$file".GQ20.hmp | cat - "$file"10K_nohdr.GQ20.hmp > "$file"10K.GQ20.hmp

# Convert to STRUCTURE format
# use SNPutilities interactively

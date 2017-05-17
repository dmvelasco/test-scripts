#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-snpid.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-snpid.txt
#SBATCH -J SNPID
#SBATCH -p bigmeml
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u


file="amygdalus"	# prefix of file name

# parse a VCF file without SNP ID and add; without ID vcf2hmp.pl script does not work properly
# Prunus VCF files using peach v1.0 genome, including scaffolds, remove word scaffold and add SNP ID

# VCF columns of interest
# position:	col1 = $1		col2 = $2		col3 = $3
# contents:	scaffold_#	position	ID (all .)
# desired:	#		position	s#_position

# perl
# remove word scaffold in file with perl one-liner (check)

# awk
# designate file separator as tab
# if the first field begins with a '#' print the line
# otherwise if the third field is a '.' then make it empty
# create a variable holding the concatenation of 's', the first field (chromosome/scaffold number),
# '_' and the second field (position) then replace the original '.' in field 3 with
# the concatenated value and print the line

zcat "$file".flt.vcf.bzip | perl -plne 's/scaffold_(\w+)/$1/' - | awk '{ OFS="\t"; if ($1 ~ /^[#]/) {print} else if ($3 ~ /\./) { sub(/\$3/, ""); temp = ("s"$1"_"$2); $3 = temp; print } ; }' - > "$file".fltID.vcf

# add -F '\t' <- only indicates input file is tab separated
# OFS ="\t" <- to designate that output file is tabl separated

# source of the complex action for the else if statement is:
# http://asdirkseesit.blogspot.com/2012/09/awk-substitute-and-concatenate-columns.html

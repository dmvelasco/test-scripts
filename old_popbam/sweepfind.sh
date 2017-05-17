#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-sweep-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-sweep-stderr.txt
#SBATCH -J sweepfind
#SBATCH -p serial
#SBATCH -a 1-8
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Declare directories
bin="/home/dmvelasc/bin"                               # program directory

i=$SLURM_ARRAY_TASK_ID

# extract information from file and add a header for sweepfinder
# ouput from PopBAM
# col 1 = chomosome
# col 2 = position
# suceeding column pairs
# frequency of SNPs at position
# number of individuals in population
# persica and dulcis are in columns 7 & 8 and 15 & 16 - determine which corresponds to which


#for i in {1..8}; do

#awk '{ OFS="\t"; print $2,$7,$8; }' all.prunus.MQ30_scaffold"$i".snp > dulcis.scaffold"$i"-temp.snp
awk '{ OFS="\t"; print $2,$17,$18; }' all.prunus.MQ30_scaffold"$i".snp > persica.scaffold"$i"-temp.snp

#echo -e "position\tx\tn" | cat - dulcis.scaffold"$i"-temp.snp > dulcis.MQ30.scaffold"$i".snp
echo -e "position\tx\tn" | cat - persica.scaffold"$i"-temp.snp > persica.MQ30.scaffold"$i".snp

#rm dulcis.scaffold"$i"-temp.snp persica.scaffold"$i"-temp.snp
rm persica.scaffold"$i"-temp.snp
#done

#for i in {1..8}; do
#"$bin"/SweepFinder -s 1000 dulcis.MQ30.scaffold"$i".snp dulcis.MQ30.scaffold"$i"_1k.sweep
#"$bin"/SweepFinder -s 1000 persica.MQ30.scaffold"$i".snp persica.MQ30.scaffold"$i"_1k.sweep

#"$bin"/SweepFinder -s 10000 dulcis.MQ30.scaffold"$i".snp dulcis.MQ30.scaffold"$i"_10k.sweep
"$bin"/SweepFinder -s 10000 persica.MQ30.scaffold"$i".snp persica.MQ30.scaffold"$i"_10k.sweep
#done

#add column with chromosome/scaffold number
# designate file separator as tab
# if the first field begins with a letter then print chrom and remainder of the line
# otherwise if the first field begins with is a number then print 's' and $i and remainder of the line
# using $i may cause a problem in the awk script because it is wrapped in single quotes

#awk '{ OFS="\t"; if ($1 ~ /^[a-z]/) { print "chrom",$1,$2,$3 } else if ($1 ~ /^[0-9]/) { print "s$i",$1,$2,$3 } ; }' dulcis.MQ30.scaffold"$i"_10k.sweep > dulcis.MQ30.s"$i"_10k.sweep
awk '{ OFS="\t"; if ($1 ~ /^[a-z]/) { print "chrom",$1,$2,$3 } else if ($1 ~ /^[0-9]/) { print "s$i",$1,$2,$3 } ; }' persica.MQ30.scaffold"$i"_10k.sweep > persica.MQ30.s"$i"_10k.sweep

# sample of sweepfinder output
# location        LR      alpha
# 38315.000000    0.352601        2.535977e-04
# 85172.769770    0.349114        2.174155e-04
# 132030.539540   0.049071        1.158977e-01
# 178888.309309   0.415672        6.015284e-04
# 225746.079079   0.018090        1.942036e-01

#GRIDSIZE
#From https://b5046268f0dde4ff9cbd60287429cd2a12b662c9.googledrive.com/host/0B9TyPcxdmHSKN3dxT19Hc1d0WUk/Programs/SweepFinder_and_Sweed.html
#The number is the number of windows, so a larger grid means a smaller window size.
#I think a window of 1 - 10 kb should be good.

#care should be taken that the grid is sufficiently dense that no large peaks appear
#when the gridsize is increased

#The significance of the maximum CLR in a region should be determined by
#analyzing data from neutral simulations which have the same sample size
#and SNP density.  When analyzing the neutral simulations, the same gridsize
#should be used as when the real data is analyzed.

#Main usage
#---------------------------------------------------------------
#The main way to use the program is with the command:
#./SweepFinder -s GRIDSIZE snp_filename out_filename
#Then, the program will read in the snps, optimize the frequency spectrum,
#and then analyze the likelihood function along a grid of size GRIDSIZE.
#The output file will have GRIDSIZE lines with the columns:
#1. location
#2. maximum CLR at location
#3. alpha corresponding to maximum CLR



#Input format
#----------------------------------------------------------------
#The snp file should be a tab-delimited file with column headers, and
#one row per SNP.  One column header should be "x" (the frequency of the
#SNP), another should be "n" (the sample size, must be greater than x), and
#another should be "position" (the chromosomal location of the SNP).
#Optionally, a fourth column named "folded" can be added.  If it is present,
#than a value of one indicates that the SNP is folded (there is no distinction
#between ancestral/derived states), and 0 means unfolded.  If the folded
#column is not present, all SNPs are assumed to be unfolded.
#Column names do not actually contain quotes.  A sample input file might look
#something like:

#position        x       n	folded
#37.000000       10      46	0
#145.000000      3       47	0
#277.000000      1       47	1
#385.000000      37      43	1
#469.000000      2       45	0
#585.000000      1       44	0
#733.000000      10      45	0


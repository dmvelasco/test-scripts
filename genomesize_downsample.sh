#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genomesize
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-genomesize.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-genomesize.txt
#SBATCH -J jelly
#SBATCH -p bigmemm
#SBATCH -a 10,20,30,40,50,60,70%2
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# number of threads
threads="8"

# Load modules
module load zlib

# Declare directories and file prefix
dir1="/home/dmvelasc/Software/jellyfish-2.1.4/bin"  	# jellyfish software binary directory
dir2="/home/dmvelasc/Software/estimate_genome_size"     # estimate genome size script directory
dir3="/group/jrigrp3/Velasco/Prunus/fastq"		# input/output directory

##### CREATE SAMPLE PREFIX #####
acc="PP01"

# make scratch directory for job
mkdir -p /scratch/dmvelasc

kmer="25"

###### JELLYFISH COUNT AND HISTOGRAM OUTPUT
# count k-mers (see jellyfish documentation for options)

# from http://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html
#jellyfish count -t 8 -C -m 25 -s 5G -o spec1_25mer --min-quality=20 --quality-start=33 */*.qf.fastq


# TWO PASS method (slower; not tested as of 2017-01-12)
#gunzip -c "$dir3"/"$acc"_1.fq.gz "$dir3"/"$acc"_2.fq.gz | "$dir1"/jellyfish bc -m "$kmer" -s 100G -t "$threads" -o "$acc"_"$kmer".bc /dev/fd/0
#gunzip -c "$dir3"/"$acc"_1.fq.gz "$dir3"/"$acc"_2.fq.gz | "$dir1"/jellyfish count -m "$kmer" -C -s 3G -t "$threads" --bc "$acc"_"$kmer".bc -o "$acc"_"$kmer".counts /dev/fd/0

# ONE PASS method
#srun gunzip -c "$dir3"/"$reads"_1.fq.gz "$dir3"/"$reads"_2.fq.gz | "$dir1"/jellyfish count -m "$kmer" -C -s 10G -t "$threads" -o "$acc"_"$kmer".counts /dev/fd/0
srun gunzip -c "$dir3"/"$acc"_1_filt_sub"$x".fq.gz "$dir3"/"$acc"_2_filt_sub"$x".fq.gz | "$dir1"/jellyfish count -m "$kmer" -C -s 6G -t "$threads" -o "$acc"_"$kmer"_sub"$x".counts /dev/fd/0

# use original reads; perl scripts later check to see that the read lengths are the same

# OTHER NOTES
# zcat does not work properly with mac os x, following is more tranportable
# gzip -dc => decompress (d) to stdout (c) OR gunzip -c is same


# Usage: jellyfish count [options] file:path+
#
# options
# -m - length of k-mer, mer length
# -o - outfile, change from default of mer_counts.jf
# -C - canonical mers of length specified
# -s - hash size (number of elements, can use letters such as M for mega- and G for giga- k-mers, i.e. 100M = 100 million elements)
# -U - skips high frequency k-mers (what is number that follows?)
# -t - number of threads to use
# /dev/fd/0 -> indicates the stdin for piping, example zcat file.fastq.gz | jellyfish count [OPTIONS] jellyfish documentation

# generate a histogram file
"$dir1"/jellyfish histo "$acc"_"$kmer"_sub"$x".counts > "$acc"_"$kmer"_sub"$x".counts.histo

# jellyfish command histo
# Usage: jellyfish histo file
# -h high setting, default 10000
# -f count all bins, including those with zero mer counts

rm "$acc"_"$kmer"_sub"$x".counts

###### ESTIMATE GENOME SIZE SCRIPTS

# generate a pdf graph of the histogram
#jellyplot.pl "$acc".counts.histo

# look at fastq.counts_0.histo.pdf and identify the approximate peak

# use find_valleys.pl to help pinpoint the actual peak
#find_valleys.pl "$acc".counts.histo

# estimate the size and coverage
#estimate_genome_size.pl --kmer=28 --peak=42 --fastq="$reads"_1.fq.gz "$reads"_2.fq.gz

# --peak - assume this value comes from jellyplot.pl and/or find_valleys.pl

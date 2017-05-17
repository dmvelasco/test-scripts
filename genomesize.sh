#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genomesize
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-genomesize.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-genomesize.txt
#SBATCH -J jelly
#SBATCH -p serial
#SBATCH -a 1-20%5
#SBATCH -t 10-00:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=24000
set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Load modules
module load zlib

# Declare directories and file prefix
dir1="/home/dmvelasc/Software/jellyfish-2.1.4/bin"  # jellyfish software binary directory
dir2="/home/dmvelasc/Software/estimate_genome_size"     # estimate genome size script directory
dir3="/home/dmvelasc/Projects/Prunus/Data/fastq"	# input/output directory

# Declare prefix array
#declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210)
#declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
# Public sequences for later - SOME HAVE MULTIPLE SRA RUNS
declare -a main=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502999 SRR502985 SRR502994 SRR502992 SRR502993 SRR502990 SRR502991 SRR502987 SRR502989 SRR502986 SRR503000 SRR503001 SRR502983 SRR502997 SRR502995 SRR502996 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
declare -a abbr=(PD11 PD12 PD13 PD14 PG01.1 PG01.2 PP01 PP02 PP03.1 PP03.2 PP04.1 PP04.2 PP05.1 PP05.2 PP06 PP07.1 PP07.2 PP08 PP09 PP10.1 PP10.2 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES
declare -a abbr=(PV03 PV04 PV05 PV06 PG02 PG04 PG05 PS04 PM01 PM02 PM03 PM04 PM05 PM06 PD11 PP37 PP39 PP40 PD16 PD17 PD18 PD21 PG03 PS03 PP38 PD20 PD19)

reads="${main["$i"]}"
acc="${abbr["$i"]}"
kmer="28"

###### JELLYFISH COUNT AND HISTOGRAM OUTPUT
# count k-mers (see jellyfish documentation for options)

# from http://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html
#jellyfish count -t 8 -C -m 25 -s 5G -o spec1_25mer --min-quality=20 --quality-start=33 */*.qf.fastq


# TWO PASS method (slower; not tested as of 2017-01-12)
#gunzip -c "$dir3"/"$reads"_1.fq.gz "$dir3"/"$reads"_2.fq.gz | "$dir1"/jellyfish bc -m "$kmer" -s 100G -t 16 -o "$acc"_"$kmer".bc /dev/fd/0
#gunzip -c "$dir3"/"$reads"_1.fq.gz "$dir3"/"$reads"_2.fq.gz | "$dir1"/jellyfish count -m "$kmer" -C -s 3G -t 16 --bc "$acc"_"$kmer".bc -o "$acc"_"$kmer".counts /dev/fd/0

# ONE PASS method
#srun gunzip -c "$dir3"/"$reads"_1.fq.gz "$dir3"/"$reads"_2.fq.gz | "$dir1"/jellyfish count -m "$kmer" -C -s 3G -t 12 -o "$acc"_"$kmer".counts /dev/fd/0
srun gunzip -c "$dir3"/"$acc"_1_filt.fq.gz "$dir3"/"$acc"_2_filt.fq.gz | "$dir1"/jellyfish count -m "$kmer" -C -s 3G -t 12 -o "$acc"_"$kmer".counts /dev/fd/0

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
# /dev/fd/0 ? is this indicating the stdin for piping, example zcat file.fastq.gz | jellyfish count [OPTIONS] jellyfish documentation

# generate a histogram file
"$dir1"/jellyfish histo "$acc"_"$kmer".counts > "$acc"_"$kmer".counts.histo

# jellyfish command histo
# Usage: jellyfish histo file
# -h high setting, default 10000
# -f count all bins, including those with zero mer counts

rm "$acc"_"$kmer".counts

###### ESTIMATE GENOME SIZE SCRIPTS

# generate a pdf graph of the histogram
#jellyplot.pl "$acc".counts.histo

# look at fastq.counts_0.histo.pdf and identify the approximate peak

# use find_valleys.pl to help pinpoint the actual peak
#find_valleys.pl "$acc".counts.histo

# estimate the size and coverage
#estimate_genome_size.pl --kmer=28 --peak=42 --fastq="$reads"_1.fq.gz "$reads"_2.fq.gz

# --peak - assume this value comes from jellyplot.pl and/or find_valleys.pl

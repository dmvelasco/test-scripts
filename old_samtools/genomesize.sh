#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Almond_BGI
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-genomesize.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-genomesize.txt
#SBATCH -J jelly
#SBATCH -p bigmemm
#SBATCH -a 1-4
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare directories and file prefix
dir1="/home/dmvelasc/bin/bin"				# jellyfish software binary directory
dir2="/home/dmvelasc/Software/estimate_genome_size"	# estimate genome size script directory

# Declare prefix array
declare -a id=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 FPS-Lovell UCD-fenzliana UCD-TNP USDA-arabica)
#PB01 PC01 PD01 PD02 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10 PD11 PD12 PD13 PD14 PF01 PG01 PK01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PP15 PR01 PS01 PS02 PT01 PU01 PV01 PV02)

declare -a fq=(FCD2GPRACXX-SZAIPI032084-64 FCC2GP0ACXX-SZAIPI032110-94_L8 FCD2H8PACXX-SZAIPI034950-14 FCD2GPRACXX-SZAIPI032087-79 FCD2H8PACXX-SZAIPI034951-15 FCD2GPRACXX-SZAIPI032085-66 FCD2GPRACXX-SZAIPI032086-74 FCD2H8PACXX-SZAIPI034949-13 FCD2GPRACXX-SZAIPI032088-80 FCD2GPRACXX-SZAIPI032089-81 FCC2GP0ACXX-SZAIPI032108-96 FCC2GP0ACXX-SZAIPI032109-95)

acc="${id["$i"]}"
reads="${fq["$i"]}"

###### JELLYFISH COUNT AND HISTOGRAM OUTPUT
# count k-mers (see jellyfish documentation for options)
gzip -dc "$acc"/"$reads"_1.fq.gz "$acc"/"$reads"_2.fq.gz | "$dir1"/jellyfish count -m 28 -o "$acc"/"$acc".counts -C -s 3G -t 16 /dev/fd/0

# use original reads; perl scripts later check to see that the read lengths are the same

# gzip -dc => decompress (d) to stdout (c)

# Usage: jellyfish count [options] file:path+
#
# options
# -m - length of k-mer, mer length
# -o - outfile, change from default of mer_counts.jf
# -C - canonical mers of length specified
# -s - hash size (number of elements, can use letters such as M for mega- and G for giga- k-mers, i.e. 100M = 100 million elements)
# -U - skips high frequecy k-mers (what is number that follows?)
# -t - number of threads to use
# /dev/fd/0 ? is this indicating the stdout <- think so, specified in jellyfish documentation

# generate a histogram
"$dir1/"jellyfish histo "$acc"/"$acc".counts > "$acc"/"$acc".counts.histo

# jellyfish subcommand of histo
# Usage: jellyfish histo file


###### ESTIMATE GENOME SIZE SCRIPTS

# generate a pdf graph of the histogram
"$dir2"/jellyplot.pl "$acc"/"$acc".counts.histo

# look at fastq.counts_0.histo.pdf and identify the approximate peak

# use find_valleys.pl to help pinpoint the actual peak
"$dir2"/find_valleys.pl "$acc"/"$acc".counts.histo

# estimate the size and coverage
#"$dir2"/estimate_genome_size.pl --kmer=28 --peak=42 --fastq="$reads"_1.fq.gz "$reads"_2.fq.gz

# --peak - assume this value comes from jellyplot.pl and/or find_valleys.pl

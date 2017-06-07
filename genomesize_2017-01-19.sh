#!/usr/bin/env bash

# Declare directories and file prefix
dir1="/usr/bin/"                           			# jellyfish software binary directory
dir2="/Volumes/SeagateX_20160628/Prunus/genomesize"     	# estimate genome size output directory
dir3="/Volumes/SeagateX_20160628/Prunus/SRAprep_postpone"	# fastq directory for some samples

# Declare prefix array
declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU1207.2 DPRU2331.9 DPRU0210)
declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
#declare -a sra=()
#declare -a pub=(PD11 PD12 PD13 PD14 PF01 PG01 PK01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PP15 PS01 PV01)

i=2

reads="${main["$i"]}"
acc="${abbr["$i"]}"

###### JELLYFISH COUNT AND HISTOGRAM OUTPUT
# count k-mers (see jellyfish documentation for options)
# TWO PASS method (slower; not tested as of 2017-01-12)
#zcat "$reads"_1.fq.gz "$reads"_2.fq.gz | "$dir1"/jellyfish bc -m 28 -s 100G -t 8 -o "$acc".bc /dev/fd/0
#zcat "$reads"_1.fq.gz "$reads"_2.fq.gz | "$dir1"/jellyfish count -m 28 -C -s 3G -t 8 --bc "$acc".bc -o "$acc".counts /dev/fd/0

# ONE PASS method
#zcat "$dir3"/"$reads"_1.fq.gz "$dir3"/"$reads"_2.fq.gz | jellyfish count -m 28 -C -s 3G -t 8 -o "$acc".counts /dev/fd/0
#zcat does not work properly with mac os x
gunzip -c "$dir3"/"$reads"_1.fq.gz "$dir3"/"$reads"_2.fq.gz | jellyfish count -m 28 -C -s 3G -t 8 -o "$acc".counts /dev/fd/0

# use original reads; perl scripts later check to see that the read lengths are the same

# gzip -dc => decompress (d) to stdout (c)

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

# generate a histogram
jellyfish histo "$acc".counts > "$acc".counts.histo

# jellyfish subcommand of histo
# Usage: jellyfish histo file
# -h high setting, default 10000
# -f count all bins, including those with zero mer counts


###### ESTIMATE GENOME SIZE SCRIPTS

# generate a pdf graph of the histogram
#jellyplot.pl "$acc".counts.histo

# look at fastq.counts_0.histo.pdf and identify the approximate peak

# use find_valleys.pl to help pinpoint the actual peak
#find_valleys.pl "$acc".counts.histo

# estimate the size and coverage
#estimate_genome_size.pl --kmer=28 --peak=42 --fastq="$reads"_1.fq.gz "$reads"_2.fq.gz

# --peak - assume this value comes from jellyplot.pl and/or find_valleys.pl

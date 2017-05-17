#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stats-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stats-stderr.txt
#SBATCH -J bamstats
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Declare program directory
dir1="/home/dmvelasc/bin"

# Load zlib 1.2.8
module load zlib

# create output files to append
#touch stats/depth.txt
#touch stats/depth_BQ20.txt stats/depth_BQ20MQ30.txt stats/depth_BQ20MQ20.txt

# flag and depth information

for i in *.bam; do
	file="${i##*/}" #gathers the file name with extension
	id="${file%.*}" #gathers only the file prefix
#	"$dir1"/samtools flagstat "$i" > stats/"$id"_stat.txt
#	"$dir1"/samtools view -bhu -q20 "$i" | "$dir1"/samtools flagstat - > stats/"$id"_stat_MQ20.txt
	awk -v x="$id" 'BEGIN{ OFS="\t"; } { print x,$0 >> "stat_MQ20.txt" NR }' stats/"$id"_stat_MQ20.txt
#	"$dir1"/samtools view -bhu -q30 "$i" | "$dir1"/samtools flagstat - > stats/"$id"_stat_MQ30.txt
#	awk -v x="$id" 'BEGIN{ OFS="\t"; } { print x,$0 >> "stat_MQ30.txt" NR }' stats/"$id"_stat_MQ30.txt
#	"$dir1"/samtools depth "$i" > stats/"$id"_depth.txt
#	awk -v x="$id" '{ sum += $3 } END { if (NR > 0) { print x"\t"sum/NR } }' stats/"$id"_depth.txt >> stats/depth.txt
#	rm stats/"$id"_depth.txt
#	"$dir1"/samtools depth -q 20 "$i" > stats/"$id"_depth.txt
#	awk -v x="$id" '{ sum += $3 } END { if (NR > 0) { print x"\t"sum/NR } }' stats/"$id"_depth.txt >> stats/depth_BQ20.txt
#	rm stats/"$id"_depth.txt
#	"$dir1"/samtools depth -q 20 -Q 20 "$i" > stats/"$id"_depth.txt
#	awk -v x="$id" '{ sum += $3 } END { if (NR > 0) { print x"\t"sum/NR } }' stats/"$id"_depth.txt >> stats/depth_BQ20MQ20.txt
#	rm stats/"$id"_depth.txt
#	"$dir1"/samtools depth -q 20 -Q 30 "$i" > stats/"$id"_depth.txt
#	awk -v x="$id" '{ sum += $3 } END { if (NR > 0) { print x"\t"sum/NR } }' stats/"$id"_depth.txt >> stats/depth_BQ20MQ30.txt
#	rm stats/"$id"_depth.txt
done

# Usage: samtools flagstat <in.bam>


# Usage: samtools depth [options] in1.bam [in2.bam [...]]
# Options:
#   -b <bed>            list of positions or regions
#   -f <list>           list of input BAM filenames, one per line [null]
#   -l <int>            minQLen
#   -q <int>            base quality threshold
#   -Q <int>            mapping quality threshold
#   -r <chr:from-to>    region

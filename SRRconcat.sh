#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fastq/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-SRRconcat.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-SRRconcat.txt
#SBATCH -J srrcat
#SBATCH -p bigmemm
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH -c 1
set -e
set -u

# Load modules
module load zlib

# 385537652 May 31 11:14 SRR502987_1.fastq.gz
#1488053746 May 31 11:15 SRR502989_1.fastq.gz
cat SRR502987_1.fastq.gz SRR502989_1.fastq.gz > SRR502987-9_1.fastq.gz
mv SRR502987-9_1.fastq.gz SRR502987_1.fastq.gz
rm SRR502989_1.fastq.gz

# 391524356 May 31 11:14 SRR502987_2.fastq.gz
#1502915271 May 31 11:15 SRR502989_2.fastq.gz
cat SRR502987_2.fastq.gz SRR502989_2.fastq.gz > SRR502987-9_2.fastq.gz
mv SRR502987-9_2.fastq.gz SRR502987_2.fastq.gz
rm SRR502989_2.fastq.gz

# 368948385 May 31 11:16 SRR502990_1.fastq.gz
#1520891721 May 31 11:16 SRR502991_1.fastq.gz
cat SRR502990_1.fastq.gz SRR502991_1.fastq.gz > SRR502990-1_1.fastq.gz
mv SRR502990-1_1.fastq.gz SRR502990_1.fastq.gz
rm SRR502991_1.fastq.gz

# 387075579 May 31 11:16 SRR502990_2.fastq.gz
#1585234261 May 31 11:17 SRR502991_2.fastq.gz
cat SRR502990_2.fastq.gz SRR502991_2.fastq.gz > SRR502990-1_2.fastq.gz
mv SRR502990-1_2.fastq.gz SRR502990_2.fastq.gz
rm SRR502991_2.fastq.gz

# 582079352 May 31 11:17 SRR502992_1.fastq.gz
#1018136567 May 31 11:18 SRR502993_1.fastq.gz
cat SRR502992_1.fastq.gz SRR502993_1.fastq.gz > SRR502992-3_1.fastq.gz
mv SRR502992-3_1.fastq.gz SRR502992_1.fastq.gz
rm SRR502993_1.fastq.gz

# 565600361 May 31 11:17 SRR502992_2.fastq.gz
#1047541359 May 31 11:18 SRR502993_2.fastq.gz
cat SRR502992_2.fastq.gz SRR502993_2.fastq.gz > SRR502992-3_2.fastq.gz
mv SRR502992-3_2.fastq.gz SRR502992_2.fastq.gz
rm SRR502993_2.fastq.gz

# 4419389035 May 31 11:22 SRR502995_1.fastq.gz
# 3450494471 May 31 11:25 SRR502996_1.fastq.gz
cat SRR502995_1.fastq.gz SRR502996_1.fastq.gz > SRR502995-6_1.fastq.gz
mv SRR502995-6_1.fastq.gz SRR502995_1.fastq.gz
rm SRR502996_1.fastq.gz

# 4073864950 May 31 11:23 SRR502995_2.fastq.gz
# 3541075333 May 31 11:26 SRR502996_2.fastq.gz
cat SRR502995_2.fastq.gz SRR502996_2.fastq.gz > SRR502995-6_2.fastq.gz
mv SRR502995-6_2.fastq.gz SRR502995_2.fastq.gz
rm SRR502996_2.fastq.gz

# 375631635 May 31 11:29 SRR502998_1.fastq.gz
#3558473835 May 31 11:31 SRR502999_1.fastq.gz
cat SRR502998_1.fastq.gz SRR502999_1.fastq.gz > SRR502998-9_1.fastq.gz
mv SRR502998-9_1.fastq.gz SRR502998_1.fastq.gz
rm SRR502999_1.fastq.gz

# 334570118 May 31 11:29 SRR502998_2.fastq.gz
#3261142700 May 31 11:32 SRR502999_2.fastq.gz
cat SRR502998_2.fastq.gz SRR502999_2.fastq.gz > SRR502998-9_2.fastq.gz
mv SRR502998-9_2.fastq.gz SRR502998_2.fastq.gz
rm SRR502999_2.fastq.gz

# 467794670 May 31 11:33 SRR503000_1.fastq.gz
#2860768713 May 31 11:34 SRR503001_1.fastq.gz
cat SRR503000_1.fastq.gz SRR503001_1.fastq.gz > SRR503000-1_1.fastq.gz
mv SRR503000-1_1.fastq.gz SRR503000_1.fastq.gz
rm SRR503001_1.fastq.gz

# 407120762 May 31 11:33 SRR503000_2.fastq.gz
#2649100980 May 31 11:35 SRR503001_2.fastq.gz
cat SRR503000_2.fastq.gz SRR503001_2.fastq.gz > SRR503000-1_2.fastq.gz
mv SRR503000-1_2.fastq.gz SRR503000_2.fastq.gz
rm SRR503001_2.fastq.gz

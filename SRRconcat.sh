#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/fastq
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-SRRconcat.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-SRRconcat.txt
#SBATCH -J srrcat
#SBATCH -p bigmemh
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
set -e
set -u

# Load modules
module load zlib

# 400776247 Apr  5 15:00 SRR502987_1.fastq.gz
#1547670376 Apr  5 15:01 SRR502989_1.fastq.gz
cat SRR502987_1.fastq.gz SRR502989_1.fastq.gz > SRR502987-9_1.fastq.gz
# 406998745 Apr  5 15:00 SRR502987_2.fastq.gz
#1564341530 Apr  5 15:02 SRR502989_2.fastq.gz
cat SRR502987_2.fastq.gz SRR502989_2.fastq.gz > SRR502987-9_2.fastq.gz

# 383631396 Apr  5 15:02 SRR502990_1.fastq.gz
#1582886250 Apr  5 15:03 SRR502991_1.fastq.gz
cat SRR502990_1.fastq.gz SRR502991_1.fastq.gz > SRR502990-1_1.fastq.gz
# 401927893 Apr  5 15:02 SRR502990_2.fastq.gz
#1649588994 Apr  5 15:04 SRR502991_2.fastq.gz
cat SRR502990_2.fastq.gz SRR502991_2.fastq.gz > SRR502990-1_2.fastq.gz

# 603547053 Apr  5 15:04 SRR502992_1.fastq.gz
#1058707601 Apr  5 15:05 SRR502993_1.fastq.gz
cat SRR502992_1.fastq.gz SRR502993_1.fastq.gz > SRR502992-3_1.fastq.gz
# 587625954 Apr  5 15:05 SRR502992_2.fastq.gz
#1090101664 Apr  5 15:06 SRR502993_2.fastq.gz
cat SRR502992_2.fastq.gz SRR502993_2.fastq.gz > SRR502992-3_2.fastq.gz

# 4678761742 Apr  5 15:11 SRR502995_1.fastq.gz
# 3610597660 Apr  5 15:15 SRR502996_1.fastq.gz
#18784285678 Apr 18 10:32 SRR502995_1.fastq.gz <- this is what it was after using
#zcat SRR502996_1.fastq.gz >> SRR502995_1.fastq.gz
cat SRR502995_1.fastq.gz SRR502996_1.fastq.gz > SRR502995-6_1.fastq.gz
# 4321728244 Apr  5 15:13 SRR502995_2.fastq.gz
# 3701167094 Apr  5 15:17 SRR502996_2.fastq.gz
#18427252180 Apr 18 10:39 SRR502995_2.fastq.gz <- this is what it was after using
#zcat SRR502996_2.fastq.gz >> SRR502995_2.fastq.gz
cat SRR502995_2.fastq.gz SRR502996_2.fastq.gz > SRR502995-6_2.fastq.gz

# 397265949 Apr  5 15:21 SRR502998_1.fastq.gz
#3768156202 Apr  5 15:23 SRR502999_1.fastq.gz
cat SRR502998_1.fastq.gz SRR502999_1.fastq.gz > SRR502998-9_1.fastq.gz
# 354868272 Apr  5 15:21 SRR502998_2.fastq.gz
#3463484651 Apr  5 15:25 SRR502999_2.fastq.gz
cat SRR502998_2.fastq.gz SRR502999_2.fastq.gz > SRR502998-9_2.fastq.gz

# 494286942 Apr  5 15:25 SRR503000_1.fastq.gz
#3027151285 Apr  5 15:27 SRR503001_1.fastq.gz
cat SRR503000_1.fastq.gz SRR503001_1.fastq.gz > SRR503000-1_1.fastq.gz
# 432338960 Apr  5 15:25 SRR503000_2.fastq.gz
#2804570910 Apr  5 15:28 SRR503001_2.fastq.gz
cat SRR503000_2.fastq.gz SRR503001_2.fastq.gz > SRR503000-1_2.fastq.gz

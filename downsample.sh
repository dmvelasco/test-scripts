#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fastq/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-seqtk.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-seqtk.txt
#SBATCH -p bigmemm
#SBATCH -J seqtk
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=8G

set -e
set -u

# Load seqtk module
module load seqtk

declare -a seed=(100 200 300 400 500 600 700)

# 10%
seqtk sample -s"${seed[0]}" PP01_1_filt.fq.gz 0.1 | gzip -c - > PP01_1_filt_sub10.fq.gz
seqtk sample -s"${seed[0]}" PP01_2_filt.fq.gz 0.1 | gzip -c - > PP01_2_filt_sub10.fq.gz

# 20%
seqtk sample -s"${seed[1]}" PP01_1_filt.fq.gz 0.2 | gzip -c - > PP01_1_filt_sub20.fq.gz
seqtk sample -s"${seed[1]}" PP01_2_filt.fq.gz 0.2 | gzip -c - > PP01_2_filt_sub20.fq.gz

# 30%
seqtk sample -s"${seed[2]}" PP01_1_filt.fq.gz 0.3 | gzip -c - > PP01_1_filt_sub30.fq.gz
seqtk sample -s"${seed[2]}" PP01_2_filt.fq.gz 0.3 | gzip -c - > PP01_2_filt_sub30.fq.gz

# 40%
seqtk sample -s"${seed[3]}" PP01_1_filt.fq.gz 0.4 | gzip -c - > PP01_1_filt_sub40.fq.gz
seqtk sample -s"${seed[3]}" PP01_2_filt.fq.gz 0.4 | gzip -c - > PP01_2_filt_sub40.fq.gz

# 50%
seqtk sample -s"${seed[4]}" PP01_1_filt.fq.gz 0.5 | gzip -c - > PP01_1_filt_sub50.fq.gz
seqtk sample -s"${seed[4]}" PP01_2_filt.fq.gz 0.5 | gzip -c - > PP01_2_filt_sub50.fq.gz

# 60%
seqtk sample -s"${seed[5]}" PP01_1_filt.fq.gz 0.6 | gzip -c - > PP01_1_filt_sub60.fq.gz
seqtk sample -s"${seed[5]}" PP01_2_filt.fq.gz 0.6 | gzip -c - > PP01_2_filt_sub60.fq.gz

# 70%
seqtk sample -s"${seed[6]}" PP01_1_filt.fq.gz 0.7 | gzip -c - > PP01_1_filt_sub70.fq.gz
seqtk sample -s"${seed[6]}" PP01_2_filt.fq.gz 0.7 | gzip -c - > PP01_2_filt_sub70.fq.gz

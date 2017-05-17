#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM2
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-ancfa-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-ancfa-stderr.txt
#SBATCH -J ancfa
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u

# Declare program directory
bin="/home/dmvelasc/bin"
dir2="/home/dmvelasc/Data/references/cerasifera"	# reference genome directory

# Load zlib 1.2.8
module load zlib

# create ancestral fasta from PC01.bam (P. cerasifera) using ANGSD -doFasta

"$bin"/angsd -i /group/jrigrp3/Velasco/Prunus/BAM2/PC01.bam -doCounts 1 -doFasta 2 -P 2 -minQ 20 -out "$dir2"/PC01
"$bin"/bwa index "$dir2"/PC01.fa.gz
gunzip -c "$dir2"/PC01.fa.gz > "$dir2"/PC01.fa
"$bin"/samtools faidx "$dir2"/PC01.fa

# angsd -doFasta
# --------------
# analysisFasta.cpp:
#	-doFasta	0
#	1: use a random base
#	2: use the most common base (needs -doCounts 1)
#	3: use the base with highest ebd (under development)
#	-minQ		13	(remove bases with qscore<minQ)
#	-basesPerLine	50	(Number of bases perline in output file)

#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-piparse.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-piparse.txt
#SBATCH -J sfsparse
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u

# %j is job allocation number


# parse pi data for P. dulcis and P. persica from popbam files
for i in {1..8}; do
	awk '{print $1,$2,$3,$4,$13,$14;}' all.prunus-10kb_"$i".sfs > dulcis-10kb_"$i".tajd.txt
        awk '{print $1,$2,$3,$4,$15,$16;}' all.prunus-10kb_"$i".sfs > dulcis-10kb_"$i".fwh.txt
        awk '{print $1,$2,$3,$4,$29,$30;}' all.prunus-10kb_"$i".sfs > persica-10kb_"$i".tajd.txt
        awk '{print $1,$2,$3,$4,$31,$32;}' all.prunus-10kb_"$i".sfs > persica-10kb_"$i".fwh.txt
done

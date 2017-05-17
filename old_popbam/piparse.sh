#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/Analysis
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-piparse.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-piparse.txt
#SBATCH -J piparse
#SBATCH -p bigmeml
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u

# %j is job allocation number


# parse pi data for P. dulcis and P. persica from popbam files
#for i in {1..8}; do
#	awk '{print $1,$2,$3,$9,$10;}' all.prunus-5kb_"$i".pi > dulcis-5kb_"$i".pi.txt
#	awk '{print $1,$2,$3,$17,$18;}' all.prunus-5kb_"$i".pi > persica-5kb_"$i".pi.txt
#        awk '{print $1,$2,$3,$9,$10;}' all.prunus-1kb_"$i".pi > dulcis-1kb_"$i".pi.txt
#        awk '{print $1,$2,$3,$17,$18;}' all.prunus-1kb_"$i".pi > persica-1kb_"$i".pi.txt
#done

# parse pi data for P. kansuensis from popbam files
for i in {1..8}; do
        awk '{print $1,$2,$3,$21,$22;}' all.prunus-5kb_"$i".pi > kansuensis-5kb_"$i".pi.txt
        awk '{print $1,$2,$3,$21,$22;}' all.prunus-1kb_"$i".pi > kansuensis-1kb_"$i".pi.txt
done

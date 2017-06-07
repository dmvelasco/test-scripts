#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fastq/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-fastqdump.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-fastqdump.txt
#SBATCH -J fastqd
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 4-00:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

# load modules
module load sratoolkit

# Create scratch directory for file download
mkdir -p /scratch/dmvelasc/

# Declare directories
fastq="/group/jrigrp3/Velasco/Prunus/fastq/"	# final directory
scratch="/scratch/dmvelasc"			# scratch directory

# other variables
# 53 SRA runs
sample="/home/dmvelasc/Projects/Prunus/Script/SRRsamples.txt"

# Declare prefix array, twenty total accession (these were the ones missed on iPlant)
declare -a main=(SRR3237762 SRR3138117 SRR3138123 SRR502984 SRR3138168 SRR3138169 SRR3237746 SRR3136181 SRR501836 SRR068359 SRR3138129 SRR3138132 SRR3138145 SRR3138147 SRR3141049 SRR3141073 SRR3141113 SRR3141181 SRR3141238 SRR3141248)

# Download SRR repositories as fastq
for x in {1..20}
do
   i=$(( x-1 ))
   reads="${main["$i"]}"
   fastq-dump -A "$reads" --split-files --defline-seq '@$sn $ri:N:0:0 length=$rl' --defline-qual '+$sn $ri:N:0:0 length=$rl' --gzip -O "$scratch"/

   # Move fastq files to final directory
   mv "$scratch"/"$reads"_1.fastq.gz "$fastq"
   mv "$scratch"/"$reads"_2.fastq.gz "$fastq"
done


while read p; do
   reads="${main["$p"]}"
   fastq-dump -A "$reads" --split-files --defline-seq '@$sn $ri:N:0:0 length=$rl' --defline-qual '+$sn $ri:N:0:0 length=$rl' --gzip -O "$scratch"/

   # Move fastq files to final directory
   mv "$scratch"/"$reads"_1.fastq.gz "$fastq"
   mv "$scratch"/"$reads"_2.fastq.gz "$fastq"
done < "$sample"

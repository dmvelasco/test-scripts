#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/selection
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-busted-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-busted-stderr.txt
#SBATCH -J busted
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of genes = number of arrays = 27864

script="/home/dmvelasc/Projects/Prunus/Script/busted_batch.sh"

#med
#c8-[62,63,65-77,86,88-96]
#c9-[65,67-70,72-77,86-97]
#c10-[64-77,86-97]
#c11-[71-77,86-97]

#bigmemm
#bigmem[1-8,10]

j=250

for i in {9111..20110.."$j"}; do
  sbatch -a "$i"-"$(( i+j ))"%1 -p med "$script"
done

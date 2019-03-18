#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-depth.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-depth.txt
#SBATCH -J depth
#SBATCH -p bigmemm
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -c 3
#SBATCH --mem=22G
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u


bin="/home/dmvelasc/bin"
BAM_dir="/group/jrigrp3/Velasco/Prunus/BAM"

# awk version may have been limited, although wasn't a problem for 2016 paper depth calculations
# but persica sample values were not the size I expected using that menthod
# bash for loop with an internal while loop took forever, only getting through 3 in 24 hours
# found below perl option; appears it may be best option
for file in "${BAM_dir}"/*sorted_markdup.bam; do
  name="${file##*/}"
  id="${name%_*.*}"
  "${bin}"/samtools depth -aa -d 0 -q 20 -Q 30 "${file}" > "${BAM_dir}"/depth_temp.txt
  awk '{print $3}' "${BAM_dir}"/depth_temp.txt > "${BAM_dir}"/temp.txt
  pos=$(cat temp.txt | wc -l)
  list="temp.txt"
  sum=`perl -nle '$sum += $_} END { print $sum' "$list"`
  avg=$(echo "scale=6 ; ${sum} / ${pos}" | bc)
  echo "${id} average coverage is ${avg}" >> markdup_BAM_depth_BQ20MQ30_dmax_perl.txt
done

#rm "${BAM_dir}"/depth_temp.txt "${BAM_dir}"/temp.txt

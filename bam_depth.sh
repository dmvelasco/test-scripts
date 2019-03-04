#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-depth.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-depth.txt
#SBATCH -J depth
#SBATCH -p bigmemh
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

#for file in "${BAM_dir}"/*HCrealign.bam; do
#  name="${file##*/}"
#  id="${name%_*.*}"
#  echo "${id}"
# all positions in genome
#  "${bin}"/samtools depth -aa "${file}" | awk -v x="${id}" '{sum += $3} END {if (NR>0) print x " average coverage is " sum/NR}' - >> "${BAM_dir}"/HCrealign_BAM_depth.txt
# all positions in genome, with base quality 20 and  mapping quality 30
#  "${bin}"/samtools depth -aa -q 20 -Q 30 "${file}" | awk -v x="${id}" '{sum += $3} END {if (NR>0) print x " average coverage is " sum/NR}' - >> "${BAM_dir}"/HCrealign_BAM_depth_BQ20MQ30.txt
# all positions in genome with base quality 20 and mapping quality 30 and no depth limitations (may be the issue with persica values)
#  "${bin}"/samtools depth -aa -d 0 -q 20 -Q 30 "${file}" | awk -v x="${id}" '{sum += $3} END {if (NR>0) print x " average coverage is " sum/NR}' - >> "${BAM_dir}"/HCrealign_BAM_depth_BQ20MQ30_dmax.txt
# NOT DONE all positions in genome (scaffolds 1-8) with base quality 20 and mapping quality 30
#  "${bin}"/samtools depth -aa -q 20 -Q 30 "${file}" | awk -v x="${id}" '{sum += $3} END {if (NR>0) print x " average coverage is " sum/NR}' - >> "${BAM_dir}"/HCrealign_BAM_depth_BQ20MQ30.txt
#done

# persica sample values still not the size I would expect, may be a limitation of awk summation
# bash for loop may be best option
for file in "${BAM_dir}"/*HCrealign.bam; do
  name="${file##*/}"
  id="${name%_*.*}"
  "${bin}"/samtools depth -aa -d 0 -q 20 -Q 30 "${file}" > "${BAM_dir}"/depth_temp.txt
  awk '{print $3}' "${BAM_dir}"/depth_temp.txt > "${BAM_dir}"/temp.txt
  pos=$(cat temp.txt | wc -l)
  list="temp.txt"
  sum=0
  while IFS='' read -r num || [[ -n "$num" ]]; do
    sum=$(( sum + num ))
  done < "$list"
#  for num in "$list"; do
#    ((sum+=num))
#  done
  avg=$(echo "scale=6 ; ${sum} / ${pos}" | bc)
  echo "${id} average coverage is ${avg}" >> HCrealign_BAM_BQ20MQ30_dmax_bashloop.txt
done

#rm "${BAM_dir}"/depth_temp.txt "${BAM_dir}"/temp.txt

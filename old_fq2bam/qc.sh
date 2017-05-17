#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Peach_GDR
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/qc-stdout-%A_%a.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/qc-stderr-%A_%a.txt
#SBATCH -a 1-23
#SBATCH -J qc
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u


module load sickle
module load zlib

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$x-1

# Declare arrays (these will need to change depending on sequence data organization)
declare -a accession=(Admiral_Dewey-DPRU1190 Babcock Bolinha Carmen-DPRU2142 Chinese_cling Diamante Dixon Dr_Davis Early_Crawford-DPRU0589 Elberta Florida_Prince-P138 Georgia_Bell JH_Hale Lovell Mayflower Nemaguard Nonpareil OHenry Okinawa Oldmixon_Free Rio_Oso_Gem Slappey-DPRU2179 St_John_Yellow-DPRU0941)
declare -a id=(PP16 PP17 PP18 PP19 PP20 PP21 PP22 PP23 PP24 PP25 PP26 PP27 PP28 PP29 PP30 PP31 PD15 PP32 PP33 PP34 PP35 PP36 PP37)

# Declare directories (also change depending on organization)
dir1="/home/dmvelasc/Software/scythe" # program directory
dir2="/group/jrigrp3/Velasco/Peach_GDR" # sequence directory prefix
dir3="/group/jrigrp3/Velasco/Prunus" # outfile directory prefix


## scythe for adapter trimming
"$dir1"/scythe -a "$dir1"/illumina_adapters.fa -o "$dir2"/scythe_"${id["$i"]}"_1.fq.gz "$dir2"/"${id["$i"]}"_1.fq.gz
"$dir1"/scythe -a "$dir1"/illumina_adapters.fa -o "$dir2"/scythe_"${id["$i"]}"_2.fq.gz "$dir2"/"${id["$i"]}"_2.fq.gz

## sickle for low quality trimming
sickle pe -f "$dir2"/scythe_"${id["$i"]}"_1.fq.gz -r "$dir2"/scythe_"${id["$i"]}"_2.fq.gz -t sanger -o "$dir2"/sickle_"${id["$i"]}"_1.fq -p "$dir2"/sickle_"${id["$i"]}"_2.fq -s "$dir2"/sickle_"${id["$i"]}"_singles.fq

gzip -c "$dir2"/sickle_"${id["$i"]}"_1.fq > "$dir2"/sickle_"${id["$i"]}"_1.fq.gz
gzip -c "$dir2"/sickle_"${id["$i"]}"_2.fq > "$dir2"/sickle_"${id["$i"]}"_2.fq.gz
gzip -c "$dir2"/sickle_"${id["$i"]}"_singles.fq > "$dir2"/sickle_"${id["$i"]}"_singles.fq.gz
rm "$dir2"/sickle_"${id["$i"]}"_1.fq "$dir2"/sickle_"${id["$i"]}"_2.fq "$dir2"/sickle_"${id["$i"]}"_singles.fq

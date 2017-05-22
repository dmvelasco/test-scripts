#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/fastq
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-GATK1prep.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-GATK1prep.txt
#SBATCH -p med
#SBATCH -J prep
#SBATCH -a 43-45
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=2600M

set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare prefix array
# Public sequences 47 total
declare -a abbr=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP12 PP14 PS01 PV01 PV03 PV04 PV05 PV06 PG02 PG04 PG05 PS04 PM01 PM02 PM03 PM04 PM05 PM06 PP37 PP39 PP40 PD16 PD17 PD18 PD21 PG03 PS03 PP38 PD20 PD19 PP11 PP13)
acc="${abbr["$i"]}"

# Create scratch directory
# -p checks to see whether it exists first; does not overwrite
mkdir -p /scratch/dmvelasc/

# Directory variables
scratch="/scratch/dmvelasc"
fq="/home/dmvelasc/Projects/Prunus/Data/fastq"

##################################################
###  Edit public fastq to have matching reads  ###
##################################################
for j in {1..2}
do
	zcat "$acc"_"$j"_filt.fq.gz | awk '{ if ($1 ~ /^@SRR/ || $1 ~ /^+SRR/) {split($1,a,"."); $1=a[1]"."a[2]; print;} else print; }' - | gzip -c - > "$scratch"/"$acc"_"$j"_filt2.fq.gz
	mv "$scratch"/"$acc"_"$j"_filt2.fq.gz "$fq"/
done

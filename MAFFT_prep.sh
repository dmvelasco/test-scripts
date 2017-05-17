#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree/fasta
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-mafft_prep.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-mafft_prep.txt
#SBATCH -J mafft
#SBATCH -a 1-27585%10
#SBATCH -p bigmemh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL

set -e
set -u

# eventual number of cycles
# -a 1-27585%50

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
#i=$(( x-1 ))

declare -a id=(PS02 PD01 PP15)
#declare -a id=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# need to modify so the array is pulled from a file with list of IDs
# also make rest of this script into a bash script that is just called by this slurm script

##### STEP 1 #####
# take the first sample in the array
# find its x-th file and extract chromosome and start position
# only needs to be done once as other files will contain the same information

sample="${id[0]}"
for match in ${x}_scaffold_*_${sample}.fa
do
	IFS='_' read -r -a array <<< "$match"
	chr="${array[2]}"
	position="${array[3]}"
	fasta="scaffold_"${chr}"_"${position}
done

##### STEP 2 #####
# make new FASTA file (name scaffold_chr_position.fa) for concatenated sequences
# determine the length of the sample array then subtract 1
# loop over all sample files for this chromosome and position
# create the sample prefix, file name, fast header name, and outfile
## Note: outfile variable can be removed if last two lines are merge
# create FASTA ID line and modify first line in FASTA
## the ID is not included in the FASTA when created by GATK from the GVCF
## appending first line in file with awk from
## http://www.unix.com/shell-programming-and-scripting/162761-append-string-first-line-file.html
# append to new FASTA file


#### NEED TO REMEMBER THAT WILL BE USING CDS FASTAS AND THEY SHOULD BE COMBINED BY GENE
#### THUS WILL NEED TO EXTRACT BY CDS, COMBINE, AND ADD GENE IDs TO THEM


touch "$fasta".fa
y="${#id[@]}"
y=$(( y-1 ))
for ((j=0; j<=y; j++)); do
	prefix="${id["$j"]}"
	filename=${x}_${fasta}_${prefix}.fa
	name=$(head -1 $filename | cut -d " " -f 2)
	awk -v pre=$prefix -v location=$name 'NR==1{$0=">"pre"_"location}1' "$filename" >> "$fasta".fa
done

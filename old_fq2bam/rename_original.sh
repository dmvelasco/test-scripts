#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/sickle
#SBATCH -o /group/jrigrp3/Velasco/Prunus/sickle/%j-stdout-vcfconvert.txt
#SBATCH -e /group/jrigrp3/Velasco/Prunus/sickle/%j-stderr-vcfconvert.txt
#SBATCH -J rename
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u


# %j is job allocation number

# renaming files from old id to new id

# Declare arrays
declare -a old=(SRR502984 SRR502998 SRR502982 DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 UCD-fenzliana USDA-arabica SRR502985 SRR502994 SRR502992 SRR502990 SRR502987 SRR502986 SRR503000 SRR502983 SRR502997 SRR502995 SRR501836 SRR068361 SRR068359 SRR068360 FPS-Lovell DPRU2578.2 UCD-TNP A_CGATGT B_TGACCA C_ACAGTG D_GCCAAT E_CAGATC F_CTTGTA G_AGTCAA H_AGTTCC)
declare -a new=(PS01 PG01 PV01 PR01 PC01 PS02 PK01 PU01 PT01 PV02 PF01 PB01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PP15 PD01 PD02 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
declare -a ext=(1 2 singles)

# Delcare directories
dir1="/group/jrigrp3/Velasco/Prunus/sickle"

x=0

for y in {0..2}

	nonsra="sickle_"${old["$x"]}"_"${ext["$y"]}".fq.gz"
	sra=""${old["$x"]}"_"${ext["$y"]}".fastq.gz"
	newfile=""${new["$i"]}"_"${ext["$y"]}".sickle.fq.gz"

	for i in $dir1; do
		if [[ $i =~ $nonsra ]]; # =~ is bash regex expression comparison operator
			then
				mv  "$nonsra" "$newfile"
				x+=

		elif [[ $i =~ $sra ]];
			then
				mv  "$sra" "$newfile"
				x+=
		else
			y=0
		fi
	done
done

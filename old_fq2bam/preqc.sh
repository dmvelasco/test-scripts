#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Peach_GDR
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/preqc-stdout-%A_%a.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/preqc-stderr-%A_%a.txt
#SBATCH -a 1-23
#SBATCH -J preqc
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u


module load zlib

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$x-1

# Declare arrays (these will need to change depending on sequence data organization)
declare -a accession=(Admiral_Dewey-DPRU1190 Babcock Bolinha Carmen-DPRU2142 Chinese_cling Diamante Dixon Dr_Davis Early_Crawford-DPRU0589 Elberta Florida_Prince-P138 Georgia_Bell JH_Hale Lovell Mayflower Nemaguard Nonpareil OHenry Okinawa Oldmixon_Free Rio_Oso_Gem Slappey-DPRU2179 St_John_Yellow-DPRU0941)
declare -a id=(PP16 PP17 PP18 PP19 PP20 PP21 PP22 PP23 PP24 PP25 PP26 PP27 PP28 PP29 PP30 PP31 PD15 PP32 PP33 PP34 PP35 PP36 PP37)

# Declare directories (also change depending on organization)
dir1="/home/dmvelasc/Projects/Prunus/Script" # program directory
dir2="/group/jrigrp3/Velasco/Peach_GDR" # sequence directory prefix
dir3="/group/jrigrp3/Velasco/Prunus" # outfile directory prefix

for a in {1..2}; do
	zcat "$dir2"/"${accession["$i"]}"/Prunus_persica-"${accession["$i"]}"."$a".fastq.gz | "$dir1"/fq2line.pl - "$dir2"/"${id["$i"]}"_"$a".fq
	# split on # and print first column with sorted read IDs then save to ID file
	awk -F "#" '{print $1 | "sort"}' "$dir2"/"${id["$i"]}"_"$a".fq > "$dir2"/"${id["$i"]}"_"$a"_readID.txt
	# if second iteration then do lots of stuff
	if [[ "$a" == 2 ]]; then
		# combine read ID files, sort, then remove unique to create master readID keys
		cat "$dir2"/"${id["$i"]}"_1_readID.txt "$dir2"/"${id["$i"]}"_2_readID.txt | sort - | uniq -d - > "$dir2"/"${id["$i"]}"_readID.txt
		cat "$dir2"/"${id["$i"]}"_1_readID.txt "$dir2"/"${id["$i"]}"_2_readID.txt | sort - | uniq - > "$dir2"/"${id["$i"]}"_readID_uniq.txt
		awk '{ print $1 "#0/1" }' "$dir2"/"${id["$i"]}"_readID.txt > "$dir2"/"${id["$i"]}"_1_readID.txt
		awk '{ print $1 "#0/2" }' "$dir2"/"${id["$i"]}"_readID.txt > "$dir2"/"${id["$i"]}"_2_readID.txt
		awk '{ print $1 "#0/1" }' "$dir2"/"${id["$i"]}"_readID_uniq.txt > "$dir2"/"${id["$i"]}"_1_readID_uniq.txt
		awk '{ print $1 "#0/2" }' "$dir2"/"${id["$i"]}"_readID_uniq.txt > "$dir2"/"${id["$i"]}"_2_readID_uniq.txt
		# remove combined read ID file
		rm "$dir2"/"${id["$i"]}"_readID.txt "$dir2"/"${id["$i"]}"_readID_uniq.txt
		# cycle through both files to find paired and unique reads
		for b in {1..2}; do
			#find paired lines and output to file: read ID file into array, check whether present in second file, output full line if true, separate with new line
			awk 'BEGIN{OFS="\n";} FNR==NR{ z[$1];next }{ for(j=1;j<=NF;j++){ if($j in z) {print $1,$2,$3,$4;} } }' "$dir2"/"${id["$i"]}"_"$b"_readID.txt "$dir2"/"${id["$i"]}"_"$b".fq | gzip -c - > "$dir2"/"${id["$i"]}"_"$b".fq.gz
			#find unpaired lines and output to file: read ID file into array, check whether not present in second file, output full line if true, separate with new line
			awk 'BEGIN{OFS="\n";} FNR==NR{ z[$1];next }{ for(j=1;j<=NF;j++){ if($j in z) {print $1,$2,$3,$4;} } }' "$dir2"/"${id["$i"]}"_"$b"_readID_uniq.txt "$dir2"/"${id["$i"]}"_"$b".fq | gzip -c - > "$dir2"/"${id["$i"]}"_"$b"_unique.fq.gz
			# replace \t with \n in read file and unique read file
			rm "$dir2"/"${id["$i"]}"_"$b".fq "$dir2"/"${id["$i"]}"_"$b"_readID.txt "$dir2"/"${id["$i"]}"_"$b"_readID_uniq.txt
		done
	else
		continue
	fi
done

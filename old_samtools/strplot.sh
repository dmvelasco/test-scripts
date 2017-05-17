#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/CLUMPP
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-structure.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-structure.txt
#SBATCH -J strplot
#SBATCH -p bigmeml
#SBATCH -a 1-6
#SBATCH -n 1
#SBATCH -c 2

set -e
set -u

# %A is array job ID
# %a is array job index

x=$SLURM_ARRAY_TASK_ID
y=$(( x+1 ))

# Declare directories
dir1="/home/dmvelasc/bin"	# program binary location
dir2="/home/dmvelasc/Projects/Prunus/Analysis/STRUCTURE/2014-10-08/harvester"
dir3="/home/dmvelasc/Projects/Prunus/Data"
# K=1 produces no CLUMPP outfile

# run in clumpp
#"$dir1"/CLUMPP -k "$y" -i "$dir2"/K"$y".indfile -o K"$y".outfile

# perl to remove colon from clumpp output file
# clumpp output file for K=2 is  something like: "1        1   (0)      1 :  0.9166 0.0834"
awk '{ OFS="\t"; s = ""; for (i = 6; i <= NF; i++) s = s $i " "; print $1"\tA\t"s }' K"$y".outfile > K"$y".out
# above from: http://stackoverflow.com/questions/5081916/how-to-print-all-the-columns-after-a-particular-number-using-awk
# modified with tab file separation

# replace genotype ID numbers with alpha-numeric genotype ID
join "$dir3"/prunus_str_id.txt K"$y".out | awk '{ s = ""; for (i = 2; i <= NF; i++) s = s $i " "; print s }' - > K"$y".txt

# format for Structure Plot
# change header, example (K=4): "group	geno	P1	P2	P3	P4"
if [ "$y" -eq 2 ]; then
	echo -e "group\tgeno\tP1\tP2" | cat - K"$y".txt > K"$y".strplot
elif [ "$y" -eq 3 ]; then
	echo -e "group\tgeno\tP1\tP2\P3" | cat - K"$y".txt > K"$y".strplot
elif [ "$y" -eq 4 ]; then
	echo -e "group\tgeno\tP1\tP2\tP3\tP4" | cat - K"$y".txt > K"$y".strplot
elif [ "$y" -eq 5 ]; then
	echo -e "group\tgeno\tP1\tP2\tP3\tP4\tP5" | cat - K"$y".txt > K"$y".strplot
elif [ "$y" -eq 6 ]; then
	echo -e "group\tgeno\tP1\tP2\tP3\tP4\tP5\tP6" | cat - K"$y".txt > K"$y".strplot
elif [ "$y" -eq 7 ]; then
	echo -e "group\tgeno\tP1\tP2\tP3\tP4\tP5\tP6\tP7" | cat - K"$y".txt > K"$y".strplot
elif [ "$y" -eq 8 ]
	echo -e "group\tgeno\tP1\tP2\tP3\tP4\tP5\tP6\tP7\tP8" | cat - K"$y".txt > K"$y".strplot
fi

# echo - e enable interpretation of backslash escapes

# header: alpha-numeric (Pop1) or alpha only (Pop)
# mandatory columns:
# first column alpha-numeric or alpha only of population IDs for each individual
# second column alpha-numeric or alpha only of individual IDs
# no extra columns
# tab delimited

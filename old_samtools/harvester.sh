#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/STRUCTURE
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-structure.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-structure.txt
#SBATCH -J structure
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 4

set -e
set -u

# %A is array job ID
# %a is array job index

# Declare directories
dir1="/home/dmvelasc/Software/structureHarvester"	# program binary location
dir2="/home/dmvelasc/Projects/Prunus/Analysis/STRUCTURE/2014-10-08/results"	# in directory
dir3="/home/dmvelasc/Projects/Prunus/Analysis/STRUCTURE/2014-10-08/harvester"	# out directory

"$dir1"/structureHarvester.py --dir="$dir2" --out="$dir3" --evanno --clump

# Usage: structureHarvester.py --dir=/path/to/results/ --out=/path/to/output
#    Options:
#      --version         show program's version number and exit
#      -h, --help        show this help message and exit
#      --dir=RESULTSDIR  The structure Results/ directory.
#      --out=OUTDIR      The out directory. If it does not exist, it will be created. Output written to summary.txt
#      --evanno          If possible, performs the Evanno 2005 method. Written to evanno.txt. default=False
#      --clumpp          Generates one K*.indfile for each value of K run, for use with CLUMPP. default=False

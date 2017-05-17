#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-map-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-map-stderr.txt
#SBATCH -J bwa2bam
#SBATCH -p serial
#SBATCH -a 1-20%5
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=24000
set -e
set -u

# %A is array job ID
# %a is array job index

##### PREREQUISITE STEPS: declare variables and arrays (check), load zlib (check), make mapping directory (check), copy reference (check)

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare prefix array
#declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210)
declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
#declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502985 SRR502994 SRR502992 SRR502990 SRR502987 SRR502986 SRR503000 SRR502983 SRR502997 SRR502983 SRR502997 SRR502995 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
#declare -a pub=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES

#reads="${main["$i"]}"
acc="${abbr["$i"]}"

# make directory to copy reference
mkdir "$acc"_temp

# Declare directories
dir1="/home/dmvelasc/bin"				# program directory
dir2="/home/dmvelasc/Data/references/bwa-peach-scf"	# reference directory
dir3="/home/dmvelasc/Projects/Prunus/Data/fastq"	# sequence directory prefix
dir4="/home/dmvelasc/Projects/Prunus/Data/BAM"		# output directory
dir5="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory

# Load zlib 1.2.8
module load zlib

# Index reference to temp directories
cp "$dir2"/Prunus_persica_v1.0_scaffolds.fa.* "$acc"_temp/
cp "$dir5"/Prunus_persica_v1.0_scaffolds.fa "$acc"_temp/

##### Map to reference with BWA MEM #####

srun "$dir1"/bwa mem -M -t 14 -k 10 "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa "$dir3"/"$acc"_1_filt.fq.gz "$dir3"/"$acc"_2_filt.fq.gz | "$dir1"/samtools view -T "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa - -o "$dir4"/"$acc".bam
# -t	threads
# -M	Mark shorter split hits as secondary (for Picard compatibility)
# -k	Minimum seed length [default 19]
# -r	Trigger re-seeding for a MEM longer than minSeedLen*FLOAT [1.5] have previously used 2.85 but may lead to lower accuracy

# remove temporary reference directory recursively
rm -r "$acc"_temp/


##### SAM AND BAM CONVERSIONS, SORTING, AND CLEAN UP #####

# sort BAM file - CONFIRM MEMORY QTY
# -@ is number of sorting and compression threads; -m is memory per thread
srun "$dir1"/samtools sort -l 9 -@ 14 -m 1G -o "$dir4"/"$acc".sorted.bam -T "$dir4"/"$acc".sorting "$dir4"/"$acc".bam

# remove unsorted BAM
rm "$dir4"/"$acc".bam

# convert back to uncompressed BAM
"$dir1"/samtools view -bh "$dir4"/"$acc".sorted.bam -o "$dir4"/"$acc".bam

# remove sorted compressed BAM, final BAM file is <prefix>.bam
rm "$dir4"/"$acc".sorted.bam


##### REMOVE DUPLICATE READS #####

# remove duplicate reads
"$dir1"/samtools rmdup -S "$dir4"/"$acc".bam "$dir4"/"$acc".nodup.bam

# remove BAM with duplicate reads
rm "$dir4"/"$acc".bam

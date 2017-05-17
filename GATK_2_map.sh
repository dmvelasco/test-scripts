#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Data/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-GATK2map-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-GATK2map-stderr.txt
#SBATCH -J map
#SBATCH -p bigmemm
#SBATCH -a 1-45%5
#SBATCH -t 16-00:00:00
#SBATCH -n 1
#SBATCH -c 4
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
#declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences
declare -a abbr=(PD11 PD12 PD13 PD14 PG01 PP01 PP02 PP03 PP04 PP05 PP06 PP07 PP08 PP09 PP10 PP11 PP12 PP13 PP14 PS01 PV01 PV03 PV04 PV05 PV06 PG02 PG04 PG05 PS04 PM01 PM02 PM03 PM04 PM05 PM06 PD11 PP37 PP39 PP40 PD16 PD17 PD18 PD21 PG03 PS03 PP38 PD20 PD19)

acc="${abbr["$i"]}"

# Declare directories
dir1="/home/dmvelasc/bin"				# program directory
dir2="/home/dmvelasc/Data/references/bwa-peach-scf"	# reference directory
dir3="/home/dmvelasc/Projects/Prunus/Data/fastq"	# sequence directory prefix
dir4="/home/dmvelasc/Projects/Prunus/Data/BAM"		# output directory
dir5="/home/dmvelasc/Data/references/persica-SCF"	# FASTA reference directory

# Load zlib 1.2.8
module load zlib

echo "create temporary directory to process sample";
date

# Make directory to copy reference
mkdir "$acc"_temp

echo "copy reference files to temporary directory";
date

# Index reference to temp directories
cp "$dir2"/Prunus_persica_v1.0_scaffolds.fa.* "$acc"_temp/
cp "$dir5"/Prunus_persica_v1.0_scaffolds.fa "$acc"_temp/


echo "Map $acc to reference with BWA MEM";
date

# Map sample to reference
# Use below with fastq that does not need the sequence ID trimmed, i.e. non-SRR samples
#srun "$dir1"/bwa mem -M -t 4 -k 10 "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa "$dir3"/"$acc"_1_filt.fq.gz "$dir3"/"$acc"_2_filt.fq.gz | "$dir1"/samtools view -T "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa - -o "$dir4"/"$acc".bam
# Use below with fastq that DOES need the sequence ID trimmed, i.e. SRR samples
srun "$dir1"/bwa mem -M -t 4 -k 10 "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa "$dir3"/"$acc"_1_filt2.fq.gz "$dir3"/"$acc"_2_filt2.fq.gz | "$dir1"/samtools view -T "$acc"_temp/Prunus_persica_v1.0_scaffolds.fa - -o "$dir4"/"$acc".bam

# -t	threads
# -M	Mark shorter split hits as secondary (for Picard compatibility)
# -k	Minimum seed length [default 19]
# -r	Trigger re-seeding for a MEM longer than minSeedLen*FLOAT [1.5]; previously used 2.85 but may lead to lower accuracy

echo "remove temp directory";
date

# remove temporary reference directory recursively
rm -r "$acc"_temp/

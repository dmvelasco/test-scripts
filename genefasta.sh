#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genomesize
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-genomesize.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-genomesize.txt
#SBATCH -J jelly
#SBATCH -p serial
#SBATCH -a 1-20%5
#SBATCH -t 10-00:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=24000
set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Load modules
module load zlib

# Declare directories and file prefix
dir1="/home/dmvelasc/Software/jellyfish-2.1.4/bin"  # jellyfish software binary directory
dir2="/home/dmvelasc/Software/estimate_genome_size"     # estimate genome size script directory
dir3="/home/dmvelasc/Projects/Prunus/Data/fastq"	# input/output directory

# Declare prefix array
#declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210)
#declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
# Public sequences for later - SOME HAVE MULTIPLE SRA RUNS
declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502999 SRR502985 SRR502994 SRR502992 SRR502993 SRR502990 SRR502991 SRR502987 SRR502989 SRR502986 SRR503000 SRR503001 SRR502983 SRR502997 SRR502995 SRR502996 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
declare -a pub=(PD11 PD12 PD13 PD14 PG01.1 PG01.2 PP01 PP02 PP03.1 PP03.2 PP04.1 PP04.2 PP05.1 PP05.2 PP06 PP07.1 PP07.2 PP08 PP09 PP10.1 PP10.2 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES

reads="${main["$i"]}"
acc="${abbr["$i"]}"

# https://www.biostars.org/p/68217/
java -Xmx2g -jar /GATK_pre/GenomeAnalysisTK-1.0.4905/GenomeAnalysisTK.jar -R reference.fasta -T FastaAlternateReferenceMaker -L sorted.bam.intervals -o consensus.fasta –variant snps_indel_file

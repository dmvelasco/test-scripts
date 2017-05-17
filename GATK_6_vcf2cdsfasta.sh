#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-vcf2fa.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-vcf2fa.txt
#SBATCH -J vcf2fa
#SBATCH -p bigmemm
#SBATCH -a 2-3
#SBATCH -t 5-00:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########

########################################################################################################
### picard verion: 2.9                                                                               ###
### GATK version: 3.7                                                                                ###
########################################################################################################

####################
### Load modules ###
####################
module load samtools/1.3.1
module load bamtools
module load java/1.8
# load GATK dependencies
module load R/3.3.1
module load maven/3.2.3
#module load GATK/3.6


#############################
### Set up the parameters ###
#############################
# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))
# full set of initial
declare -a id=(PS02 PD01 PP15)
sample="${id["$i"]}"

# location for the picard.jar
picard="/home/dmvelasc/Software/picard/picard.jar"
# location for the GenomeAnalysisTK.jar
GATK="/home/dmvelasc/Software/GATK/GenomeAnalysisTK.jar"
# genome reference file location
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# file directory for $sample vcf to convert to fasta
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK"
# gene intervals
intervals="/home/dmvelasc/Data/references/persica-SCF/cds_intervals"

# Declare directories and file prefix
dir1="/home/dmvelasc/Software/bedtools2/bin"		# bedtools location
dir2="/home/dmvelasc/Projects/Prunus/Data/BAM"		# BAM directory
dir3="/home/dmvelasc/Data/references/persica-SCF"	# gff3 location
dir4="/scratch/dmvelasc"

####################
### Begin script ###
####################
echo "begin CDS FASTA script"
date

# create scratch directory for temporary file placement
mkdir -p /scratch/dmvelasc/"$sample"

# loop to create and process each CDS

while read p; do
	# GATK portion - create CDS FASTA
	java -Xmx2g -jar "$GATK" \
	-R "$genome" \
	-T FastaAlternateReferenceMaker \
	-o "$dir4"/"$sample"/"$p"_"$sample"_cds.fa \
	-IUPAC "$sample" \
	-raw \
	--variant "$vcf"/"$sample".g.vcf \
	-L "$intervals"/"$p".intervals
	# basic manipulations to create final CDS FASTA
	echo ">${p}_${sample}" > "$dir4"/"$sample"/"$p"_"$sample".fa
	awk '{printf $0;}' "$dir4"/"$sample"/"$p"_"$sample"_cds.fa | fold -w 60 - >> "$dir4"/"$sample"/"$p"_"$sample".fa
	rm "$dir4"/"$sample"/"$p"_"$sample"_cds.fa
done < "$dir3"/Prunus_persica_v1.0_genes_list.gff3

# move sample file directory from scratch
mv /scratch/dmvelasc/"$sample"/ /home/dmvelasc/Projects/Prunus/Analysis/genetree/

echo "end CDS FASTA script"
date

#################################
### Other notes and resources ###
#################################

# from: https://www.biostars.org/p/17705/
# still how to get specific regions..., -L option with regions from gff?
# does -L take region file? sort of: -T SelectVariants \ -L /path/to/my.interval_list
# from: http://gatkforums.broadinstitute.org/gatk/discussion/2441/can-selectvariants-be-used-to-limit-vcf-files-by-interval-list
# see: https://software.broadinstitute.org/gatk/documentation/article?id=1319
# or overlap after with bedtools?

# GATK-accepted interval lists
# https://software.broadinstitute.org/gatk/documentation/article?id=1319

#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fasta/fasta-aln
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-macse2.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-macse2.txt
#SBATCH -J macse2
#SBATCH -a 41-4000%10
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of total array jobs
# -a 1-27864

# Load modules
module load java

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare directories
macse="/home/dmvelasc/Software/macse2/macse_v2.01.jar"		# macsev2 program
ref="/home/dmvelasc/Data/references/persica-SCF"		# reference directory
scratch="/scratch/dmvelasc/fasta-pep"				# scratch directory
dna_dir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-concat"	# directory of CDS fasta dna sequences
dna_aln="/group/jrigrp3/Velasco/Prunus/fasta/fasta-macse2"	# directory of CDS aligned fasta backtranslated dna sequences
dna_aln_final="/group/jrigrp3/Velasco/Prunus/fasta/fasta-nostop"
pep_dir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-peptide"	# amino acid alignment directory

####################
### Begin script ###
####################

# multi-sequence alignment with MACSE v2
# translates nucleotide sequence to amino acids and simultaneously aligns both
# outputs both alignments

#### STEP 1: CAPTURE GENE ID
mapfile -s "$i" -n 1 -t gene < "$ref"/Prunus_persica_v1.0_genes_list.gff3
echo -e "begin MACSE v2 alignment and export for CDS sequence ${gene[0]}"
date

#### STEP 2: ALIGN MULTI-SEQUENCE DNA FASTA
# sequence alignment of nucleotides and translated amino acid sequences
java -Xmx30G -jar "$macse" -prog alignSequences -seq "$dna_dir"/"${gene[0]}"_cds.fa -out_NT "$dna_aln"/"${gene[0]}"_cds_aln.fasta -out_AA "$pep_dir"/"${gene[0]}"_cds_pep_aln.fasta

#### STEP 3: EXPORT NUCLEOTIDE ALIGNMENT FOR BUSTED
# modify alignment to insert gaps at frameshifts and stop codons for use with HyPhy's BUSTED
java -Xmx30G -jar "$macse" -prog exportAlignment -align "$dna_aln"/"${gene[0]}"_cds_aln.fasta -charForRemainingFS - -codonForExternalFS --- -codonForFinalStop --- -codonForInternalFS --- -codonForInternalStop --- -out_NT "$dna_aln_final"/"${gene[0]}"_cds_aln_nostop.fasta -out_stat_per_site "$dna_aln_final"/"${gene[0]}"_cds_aln_stats.csv

echo -e "end MACSE v2 alignment and export for CDS sequence ${gene[0]}"
date

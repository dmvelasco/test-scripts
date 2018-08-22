#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fasta/fasta-concat
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-tcoffee_run.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-tcoffee_run.txt
#SBATCH -J tcoffee
#SBATCH -a 1-10%2
#SBATCH -p bigmemh
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=bigmem1
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of total array jobs
# -a 1-27864%50

# Load zlib 1.2.8
module load zlib
module load emboss
module load blast

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

# Declare directories
ref="/home/dmvelasc/Data/references/persica-SCF"		# reference directory
scratch="/scratch/dmvelasc/fasta-pep"				# scratch directory
dna_dir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-concat"	# directory of CDS fasta dna sequences
dna_aln_dir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-alntx"	# directory of CDS aligned fasta backtranslated dna sequences
pep_dir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-pep"		# directory of CDS fasta peptide sequences
pep_aln_dir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-pep"	# directory of CDS aligned fasta peptide sequences

####################
### Begin script ###
####################

# translate multi-sequence FASTA then align and back translate to DNA alignment
# using original multi-sequence DNA FASTA
# all with T-Coffee program and features

#### STEP 1: CAPTURE GENE ID
mapfile -s "$i" -n 1 -t gene < "$ref"/Prunus_persica_v1.0_genes_list.gff3
echo -e "begin T-coffee cycle for CDS sequence ${gene[0]}"
date

#### STEP 2: TRANSLATE MULTI-SEQUENCE DNA FASTA
t_coffee -other_pg seq_reformat -in "$dna_dir"/"${gene[0]}"_cds.fa -action +translate -output fasta_seq > "$pep_dir"/"${gene[0]}"_cds_pep.fasta

#### STEP 3: ALIGN MULTI-SEQUENCE PEPTIDE FASTA
t_coffee "$pep_dir"/"${gene[0]}"_cds_pep.fasta -mode expresso -blast_server=NCBI

#### STEP 4: BACK TRANSLATE MULTI-SEQUENCE PEPTIDE FASTA WITH ORIGINAL MULTI-SEQUENCE DNA FASTA
#t_coffee -other_pg seq_reformat -in "$dna_dir"/"${gene[0]}"_cds.fa -in2 "$pep_aln_dir"/"${gene[0]}"_cds_pep.aln -action +thread_dna_on_prot_aln -output fasta_aln > "$dna_aln_dir"/"${gene[0]}"_cds.aln

echo -e "end T-Coffee cycle for CDS sequence ${gene[0]}"
date

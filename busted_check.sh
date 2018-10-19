#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/selection
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-check-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-check-stderr.txt
#SBATCH -J check
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p bigmemm
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of genes = number of arrays = 27864
# --exclude=bigmem1,bigmem2

### Declare directories ###
ref="/home/dmvelasc/Data/references/persica-SCF"		# reference directory
gene_list="${ref}/Prunus_persica_v1.0_genes_list.gff3"		# gene list
seldir="/home/dmvelasc/Projects/Prunus/Analysis/selection"	# final fasta directory

for i in {0..27863}; do
  #### sample ID file
  mapfile -s "$i" -n 1 -t gene < "${gene_list}"
  # declare gene ID variable from array
  gene_id="${gene[0]}"

  j=$(( i+1))

  ##### B U S T E D #####
  if [ ! -f "${seldir}/${gene_id}_cds_aln_nostop.fasta.BUSTED.json" ]; then
  # if file does not exist
    echo "${j} : ${gene_id}_cds_aln_nostop.fasta.BUSTED.json does not exist"
  fi
done

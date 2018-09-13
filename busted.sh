#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/selection
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-busted-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-busted-stderr.txt
#SBATCH -J busted
#SBATCH -p bigmemm
#SBATCH -a 1-40%1
#SBATCH -t 8-00:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --exclude=bigmem1,bigmem2
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of genes = number of arrays = 27864

### Load modules ###
module load hyphy

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
ref="/home/dmvelasc/Data/references/persica-SCF"				# reference directory
gene_list="${ref}/Prunus_persica_v1.0_genes_list.gff3"				# gene list
FASTAdir="/group/jrigrp3/Velasco/Prunus/fasta/fasta-nostop"			# final fasta directory
script="/home/dmvelasc/Projects/Prunus/Script"					# script directory that includes batch file and consensus tree
busted="/share/apps/hyphy-2.3.13/lib/TemplateBatchFiles/SelectionAnalyses"	# BUSTED batch file directory, trying to overcome misplacement of libv3 directory
finaldir="/home/dmvelasc/Projects/Prunus/Analysis/selection"
#### sample ID file
# column 1: ID, column2: other ID/information
#list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

echo -e "extract sample ID from Script/sample.txt using mapfile"
date

mapfile -s "$i" -n 1 -t gene < "${gene_list}"
# declare gene ID variable from array
gene_id="${gene[0]}"
echo -e "Starting BUSTED analysis for ${gene_id}"
date

##### B U S T E D #####
if [ -s "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta" ]
# if file exists and is not empty
then
  #HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${busted}"/BUSTED.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta" "${script}/splitstree_all.tree" "All" "" > "$finaldir"/"${gene_id}_cds_busted.txt"
  HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${busted}"/BUSTED.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta" "${script}/split_alltest2_final_4_rooted_busted.nwk" "All" "" > "$finaldir"/"${gene_id}_cds_busted_frgnd.txt"
  if [ -s "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta.BUSTED.json" ]
  then
    mv "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta.BUSTED.json" "$finaldir"
  fi
  #HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${busted}"/BUSTED.bf "Universal" "${FASTAdir}/${gene_id}_cds.aln" "${script}/split_alltest2_final_4_rooted.nwk" "All" "" > "$finaldir"/"${gene_id}_cds_busted_none.txt"
  # what is output directory? same as input directory. what is output file? BUSTED.json
fi
done

echo -e "BUSTED analysis finished"
date

#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/selection
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-absrel-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-absrel-stderr.txt
#SBATCH -J absrel
#SBATCH -a 1-20%5
#SBATCH -t 8-00:00:00
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

# number of genes with evidence of selection from BUSTED = number of arrays = 4736
# --exclude=bigmem1,bigmem2

### Load modules ###
module load hyphy

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))

### Declare directories ###
ref="/home/dmvelasc/Data/references/persica-SCF"					# reference directory
gene_list="/home/dmvelasc/Projects/Prunus/Analysis/selection/busted_p-values.txt"	# gene list
FASTAdir="/home/dmvelasc/Projects/Prunus/Data/fasta/fasta-nostop"			# final fasta directory
script="/home/dmvelasc/Projects/Prunus/Script"						# script directory that includes batch file and consensus tree
batchdir="/share/apps/hyphy-2.3.13/lib/TemplateBatchFiles/SelectionAnalyses"		# BUSTED batch file directory
finaldir="/home/dmvelasc/Projects/Prunus/Analysis/select_branches"
#### sample ID file
# column 1: ID, column2: other ID/information
#list="/home/dmvelasc/Projects/Prunus/Script/sample.txt"

echo -e "extract sample ID from ${script}/${gene_list} using mapfile"
date

mapfile -s "$i" -n 1 -t id < "${gene_list}"
# declare gene ID variable from array
gene=(`echo "{id[0]}"`)
gene_id="${gene[0]}"
echo -e "Starting aBSREL analysis for ${gene_id}"
date

##### a B s R e L #####
if [ -s "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta" ]
# if file exists and is not empty
then
  #HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batchdir}"/aBSREL.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta" "${script}/splitstree_all.tree" "All" "" > "$finaldir"/"${gene_id}_cds_absrel.txt"
  HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batchdir}"/aBSREL.bf "Universal" "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta" "${script}/split_alltest2_final_4_rooted_busted.nwk" "All" "" "${finaldir}/${gene_id}_cds_absrel.out" > "${finaldir}/${gene_id}_cds_absrel.txt"
#  if [ -s "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta.BUSTED.json" ]
#  then
#    mv "${FASTAdir}/${gene_id}_cds_aln_nostop.fasta.BUSTED.json" "$finaldir"
#  fi
  #HYPHYMP LIBPATH=/share/apps/hyphy-2.3.13/lib/ "${batchdir}"/aBSREL.bf "Universal" "${FASTAdir}/${gene_id}_cds.aln" "${script}/split_alltest2_final_4_rooted.nwk" "All" "" > "$finaldir"/"${gene_id}_cds_absrel_none.txt"
  # what is output directory? same as input directory. what is output file? BUSTED.json
fi

echo -e "aBSREL analysis finished"
date

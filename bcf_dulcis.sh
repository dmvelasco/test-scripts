#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf.txt
#SBATCH -J vcf
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u


# %A is array job ID
# %a is array job index
# %j is job allocation number

#x=$SLURM_ARRAY_TASK_ID
#i=$(( x-1 ))

# Load zlib 1.2.8
module load zlib

# Declare BAM file list array
#declare -a group=(dulcis persica almond peach amygdalus prunus)

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Script/old_samtools"	# script directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory

##### PILEUP, CALL VARIANTS, FILTER VARIANTS #####
# PILEUP OF SORTED, UNCOMPRESSED BAM TO BCF
# filter VCF (approximate average depth is 30X)

# pileup and filtering of genomic SNPs
#"$dir1"/samtools mpileup -BRug -t DP,SP -f "$dir3"/Prunus_persica_v1.0_scaffolds.fa -b "$dir2"/all.dulcisbam.txt | "$dir1"/bcftools call -O z -Avm -f GQ,GP - > dulcis.raw.vcf.bzip
#zcat dulcis.raw.vcf.bzip | perl -plne 's/scaffold_(\w+)/$1/' - | bgzip -cf > dulcis.noscaf.vcf.bzip
#"$dir1"/bcftools index -f dulcis.raw.vcf.bzip
#"$dir1"/bcftools view -O z -o dulcis.flt_s5.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %MIN(GQ)>=30" "-r" scaffold_5:10900000-13700000 dulcis.raw.vcf.bzip
#"$dir1"/bcftools view -O z -o dulcis.flt_s5.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %MIN(GQ)>=30" "-r" scaffold_5:10000000-15000000 dulcis.raw.vcf.bzip
"$dir1"/bcftools view -O z -o dulcis.raw_s5.vcf.bzip "-r" scaffold_5:10000000-15000000 dulcis.raw.vcf.bzip

####the below line works filtering seems too stringent for above region (with %AVG(DP)<=100)
#"$dir1"/bcftools view -O z -o dulcis.flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %AVG(DP)<=100 && %MIN(GQ)>=30" "-r" 5 dulcis.noscaf.vcf.bzip

#"$dir1"/bcftools view -O z -o dulcis.flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %MIN(GQ)>=30" dulcis.raw.vcf.bzip

#"$dir1"/bcftools index -f almond.raw.vcf.bzip
#"$dir1"/bcftools view -O z -o almond.flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %MIN(GQ)>=30" almond.raw.vcf.bzip


#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf.txt
#SBATCH -J vcf
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 8
set -e
set -u


#x=$SLURM_ARRAY_TASK_ID
#i=$(( x-1 ))

# Load zlib 1.2.8
module load zlib

# Declare BAM file list array
#declare -a group=(dulcis persica almond peach amygdalus prunus)

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"

##### PILEUP, CALL VARIANTS, FILTER VARIANTS #####
# PILEUP OF SORTED, UNCOMPRESSED BAM TO BCF
# filter VCF (approximate average depth is 15-30X)

# pileup and filtering of genomic SNPs
# original, note call has -A (keeps all variants found even if not present in final genotypes) and -v (variants only)
#"$dir1"/samtools mpileup -BRug -t DP,SP -f "$dir3"/Prunus_persica_v1.0_scaffolds.fa -b all."${group["$i"]}"bam.txt | "$dir1"/bcftools call -O z -Avm -f GQ,GP - > "${group["$i"]}".raw.vcf.bzip
#"$dir1"/bcftools index "${group["$i"]}".raw.vcf.bzip
#"$dir1"/bcftools view -O z -o "${group["$i"]}".flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %AVG(DP)<=500 && %MIN(GQ)>=30" "${group["$i"]}".raw.vcf.bzip
#"$dir1"/bcftools view -O z -o "${group["$i"]}".flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>5 && %MIN(GQ)>=30" prunus.raw.vcf.bzip
#"$dir1"/bcftools view -O z -o "${group["$i"]}".flt_test.vcf.bzip -s PD03.bam,PP04.bam,PR01.bam,PU01.bam,PV02.bam,PT01.bam,PC01.bam -i "%QUAL>=30 && %AVG(DP)>5 && %MIN(GQ)>=30" "${group["$i"]}".raw.vcf.bzip
# with -A removed (consider removing -v, which only calls variant sites)
"$dir1"/samtools mpileup -BRug -t DP,SP -f "$dir3"/Prunus_persica_v1.0_scaffolds.fa -b all.prunusbam.txt | "$dir1"/bcftools call -O z -vm -f GQ,GP - > "$dir2"/prunus-A.raw.vcf.bzip
"$dir1"/bcftools index "$dir2"/prunus-A.raw.vcf.bzip
#"$dir1"/bcftools view -O z -o "$dir2"/prunus-A.flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>3 && %AVG(DP)<=500 && %MIN(GQ)>=30" "$dir2"/prunus-A.raw.vcf.bzip
"$dir1"/bcftools view -O z -o "$dir2"/prunus-A.flt.vcf.bzip -i "%QUAL>=30 && %AVG(DP)>5 && %MIN(GQ)>=30" "$dir2"/prunus.raw.vcf.bzip
"$dir1"/bcftools view -O z -o "$dir2"/prunus-A.flt_test.vcf.bzip -s PD03.bam,PP04.bam,PR01.bam,PU01.bam,PV02.bam,PT01.bam,PC01.bam -i "%QUAL>=30 && %AVG(DP)>5 && %MIN(GQ)>=30" "$dir2"/prunus-A.raw.vcf.bzip

#!/bin/bash
#SBATCH -D /group/jrigrp3/Velasco/Prunus/BAM/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-vcf.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-vcf.txt
#SBATCH -J qry
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -c 2
set -e
set -u

# Load zlib 1.2.8
module load zlib

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"

##### BCFtools options for extracting information from vcf files #####
# query

# query to extract chromosome, position, and genotype (both 0,1,2 and ACGT, IUPAC)
#"$dir1"/bcftools view -g het "$dir2"/prunus-A.flt_test.vcf.bzip | "$dir1"/bcftools query -o "$dir2"/prunus-A.query_test.vcf.gz -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%TGT\t%GQ]\n' -
#"$dir1"/bcftools query -o "$dir2"/prunus-A.query_test.vcf.gz -i "AC>0" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%TGT\t%GQ]\n' "$dir2"/prunus-A.flt_test.vcf.bzip
"$dir1"/bcftools query -o "$dir2"/prunus-A.query_testIUPAC.vcf.gz -i "AC>0" -f '%CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%IUPACGT\t%GQ]\n' "$dir2"/prunus-A.flt_test.vcf.bzip

# add header line of...
echo "CHROM      POS        TYPE	REF     ALT     PD_03        GQ      PP_04        GQ      PR_01        GQ      PU_01        GQ      PV_02        GQ      PT_01        GQ      PC_01        GQ" > "$dir2"/temp.txt
grep SNP "$dir2"/prunus-A.query_testIUPAC.vcf.gz >> "$dir2"/temp.txt
awk '{OFS="\t";} {print $6,$8,$10,$12,$14,$16,$18;}' "$dir2"/temp.txt > "$dir2"/prunus-A_SNAPP_test.txt


# Transpose columns and rows
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END { 
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' "$dir2"/prunus-A_SNAPP_test.txt > "$dir2"/prunus-A_SNAPP_test_xpose.txt

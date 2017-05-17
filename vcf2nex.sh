#!/bin/bash
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/VCF
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-snapp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-snapp.txt
#SBATCH -J snapp
#SBATCH -p bigmemm
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 0:30:00
set -e
set -u

# Load zlib 1.2.8
module load zlib
module load java

# Declare directories
dir1="/home/dmvelasc/bin"					# software binary directory
dir2="/home/dmvelasc/Projects/Prunus/Analysis/VCF"		# VCF directory
dir3="/home/dmvelasc/Data/references/persica-SCF"		# FASTA reference directory
dir4="/group/jrigrp3/Velasco/Prunus/BAM"

# Declare other variables
infile="prunus-A.flt.vcf.bzip"		# starting vcf
out1="amyg_split"			# filter on quality, depth, genotype quality
out2="amyg_split_noscaf"		# query for genotype and genotype quality
thin="5000"				# spacing between each SNP
out3="amyg_split_thin5000"		# "scaffold_" removed
out4="amyg_split_query"			# final thinned vcf
bcfquery="split_amygtest.txt"		# file from bcf query
xpose="split_amygtest_final.txt"	# transposed file
outfile="split_amygtest_final.nex"	# final file name

# NEXUS variables
#DIMENSIONS (see below)
#FORMAT
datatype="standard" #from Integer, which is what SNAPP uses
symbols='"012"'
interleave="no"
missing="?"

##### BCFtools options for extracting information from vcf files #####
# filter vcf for chromosomes 1-8, snps only
#"$dir1"/bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 -s PD03.bam,PP04.bam,PR01.bam,PU01.bam,PV02.bam,PT01.bam -i "AC>0 && %QUAL>=250 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out1".vcf -v snps "$dir2"/"$infile"
#"$dir1"/bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 -s PD03.bam,PP04.bam,PU01.bam,PV02.bam,PT01.bam -i "AC>0 && %QUAL>=250 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out1".vcf -v snps "$dir2"/"$infile"

"$dir1"/bcftools index -f "$infile"

#"$dir1"/bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 -s PD03.bam,PD04.bam,PD05.bam,PD09.bam,PD10.bam,PD06.bam,PP02.bam,PP04.bam,PP03.bam,PP05.bam,PP11.bam,PP13.bam,PU01.bam,PV01.bam,PV02.bam,PS01.bam,PS02.bam,PB01.bam,PG01.bam,PT01.bam -i "AC>0 && %QUAL>=250 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out1".vcf -v snps "$dir2"/"$infile"

# below for splitstree
"$dir1"/bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8 -i "AC>0 && %QUAL>=250 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out1".vcf -v snps "$dir2"/"$infile"

#"$dir1"/bcftools view -r 1,2,3,4,5,6,7,8 -i "AC>0 && %QUAL>=250 && %AVG(DP)>5 && %MIN(GQ)>=30" -O v -o "$out2".vcf -v snps "$dir2"/"$out1".vcf
# -i, --include EXPRESSION
# above for variant quality, depth, and genotype quality
# -r, --regions chr|chr:pos|chr:from-to|chr:from-[,…]
# Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.

# remove "scaffold_"
perl -plne 's/scaffold_(\w+)/$1/' "$out1.vcf" > "$out2".vcf

# thin to minimum spacing between SNPs
"$dir1"/vcftools --vcf "$out2".vcf --thin "$thin" --recode --out "$out3"

# query for only genotypes and genotype qualities
"$dir1"/bcftools query -o "$dir2"/"$out4".vcf -f '%CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT\t%GQ]\n' "$dir2"/"$out3".recode.vcf

# add header line of...
#echo -e \t "CHROM\tPOS\tTYPE\tREF\tALT\tPD_03 \tGQ\tPP_04 \tGQ\tPR_01 \tGQ\tPU_01 \tGQ\tPV_02 \tGQ\tPT_01 \tGQ" > "$dir2"/temp.txt
#echo -e \t "CHROM\tPOS\tTYPE\tREF\tALT\tPD_03 \tGQ\tPP_04 \tGQ\tPU_01 \tGQ\tPV_02 \tGQ\tPT_01 \tGQ" > "$dir2"/temp.txt
echo -e \t "CHROM\tPOS\tTYPE\tREF\tALT\tPD_03 \tGQ\tPD_04 \tGQ\tPD_05 \tGQ\tPD_09 \tGQ\tPD_10 \tGQ\tPD_06 \tGQ\tPP_02 \tGQ\tPP_04 \tGQ\tPP_03 \tGQ\tPP_05 \tGQ\tPP_11 \tGQ\tPP_13 \tGQ\tPU_01 \tGQ\tPV_01 \tGQ\tPV_02 \tGQ\tPS_01 \tGQ\tPS_02 \tGQ\tPB_01 \tGQ\tPG_01 \tGQ\tPT_01 \tGQ" > "$dir2"/temp.txt
#PD03,PD04,PD05,PD09,PD10,PD06,PP02,PP04,PP03,PP05,PP11,PP13,PU01,PV01,PV02,PS01,PS02,PB01,PG01,PT01

# select bialellic SNPs only by selecting lines that DO NOT include a comma
grep -v , "$out4".vcf >> "$dir2"/temp.txt

# select only genotype columns
#awk -F"\t" '{OFS="\t";} {print $6,$8,$10,$12,$14,$16,$18;}' "$dir2"/temp.txt > "$dir2"/"$bcfquery"
#awk -F"\t" '{OFS="\t";} {print $6,$8,$10,$12,$14,$16;}' "$dir2"/temp.txt > "$dir2"/"$bcfquery"
#awk -F"\t" '{OFS="\t";} {print $6,$8,$10,$12,$14;}' "$dir2"/temp.txt > "$dir2"/"$bcfquery"
awk -F"\t" '{OFS="\t";} {print $6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44;}' "$dir2"/temp.txt > "$dir2"/"$bcfquery"

# convert genotypes to 0,1,2 format with 0 homozygous reference (0/0), 1 heterozygous (0/1), 2 as homozygous alternate (1/1)
# substitute 0 for 0/0, 1 for 0/1, 2 for 1/1, and . for ./. (these are missing data, no missing data present)
# 1 at end of awk command prints line
awk -F"\t" '{OFS="\t";} {gsub(/0\/0/,"0");gsub(/0\/1/,"1");gsub(/1\/1/,"2");gsub(/.\/./,".")}1' "$dir2"/"$bcfquery" > temp_xpose.txt

# Transpose columns and rows
awk -F"\t" '
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
            str=str""a[i,j];
        }
        print str
    }
}' "$dir2"/temp_xpose.txt > "$dir2"/"$xpose"
#above original was str=str" "a[i,j];


# Nexus file variables requiring intermediate files
Ntax=$(cat "$xpose" | wc -l)	#number of taxa
Nchar=$(( `cat "$bcfquery" | wc -l`- 1 )) #number of snps

# Nexus file header and body creation
echo "#NEXUS">"$outfile"
echo "[Written $(date)]">>"$outfile"
echo "BEGIN Data;">>"$outfile"
echo -e "\tDIMENSIONS NTAX="$Ntax" NCHAR="$Nchar";">>"$outfile"
echo -e "\tFORMAT DATATYPE="$datatype" Symbols="$symbols" INTERLEAVE="$interleave" missing="$missing";" >> "$outfile"
echo "Matrix">>"$outfile"
cat "$xpose">>"$outfile"
echo ";">>"$outfile"
echo "END;">>"$outfile"

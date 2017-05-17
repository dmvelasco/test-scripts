#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/summstats
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-intersect-stdout.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-intersect-stderr.txt
#SBATCH -J intersect
#SBATCH -p bigmeml
#SBATCH -n 1
#SBATCH -c 4
set -e
set -u

# Directory variables
dir1="/home/dmvelasc/Data/references/persica-SCF" # program directory
gff3="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_genes_only_sort.gff3"

# header for out files
echo -e "chr\twinstart\twinstop\tgenestart\tgenestop\tintersect\tstrand\tthetaW\tthetaPi\tTajima\tfayh\tzeng\tgeneID" > intersect_header.txt
echo -e "chr\twinstart\twinstop\tgenestart\tgenestop\tintersect\tstrand\tunweighted\tweighted\tgeneID" > intersect_fst_header.txt

# further prepped gff3 file with following
# perl -plne 's/scaffold_(\w+)/$1/' Prunus_persica_v1.0_genes_only_sort.gff3  > test.gff3
# was creating error because query file only had numbers for chromosomes whereas gff3 file had prefix of "scaffold_"

# bedtools to find intersection of zeng windows and peach gff file
# split gene information to get only gene ID, select columns of interest, sort, and concatenate with head to final file

# Zeng's E
bedtools intersect -sorted -wao -a PD_windows_z5.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PD_intersect_z5.txt
bedtools intersect -sorted -wao -a PP_windows_z5.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PP_intersect_z5.txt

bedtools intersect -sorted -wao -a PD_windows_z01.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PD_intersect_z01.txt
bedtools intersect -sorted -wao -a PP_windows_z01.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PP_intersect_z01.txt

bedtools intersect -sorted -wao -a PD_windows_z001.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PD_intersect_z001.txt
bedtools intersect -sorted -wao -a PP_windows_z001.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PP_intersect_z001.txt

# Theta Pi
bedtools intersect -sorted -wao -a PD_windows_pi5.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PD_intersect_pi5.txt
bedtools intersect -sorted -wao -a PP_windows_pi5.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PP_intersect_pi5.txt

bedtools intersect -sorted -wao -a PD_windows_pi1.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PD_intersect_pi1.txt
bedtools intersect -sorted -wao -a PP_windows_pi1.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PP_intersect_pi1.txt

bedtools intersect -sorted -wao -a PD_windows_pi01.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PD_intersect_pi01.txt
bedtools intersect -sorted -wao -a PP_windows_pi01.txt -b "$gff3" | awk '{OFS="\t"; split($17, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$12,$13,$18,$15,$4,$5,$6,$7,$8,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_header.txt - > PP_intersect_pi01.txt

# Fst
bedtools intersect -sorted -wao -a PE_windows_fst1.txt -b "$gff3" | awk '{OFS="\t"; split($14, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$9,$10,$15,$12,$4,$5,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_fst_header.txt - > PE_intersect_fst1.txt
bedtools intersect -sorted -wao -a PE_windows_fst05.txt -b "$gff3" | awk '{OFS="\t"; split($14, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$9,$10,$15,$12,$4,$5,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_fst_header.txt - > PE_intersect_fst05.txt
bedtools intersect -sorted -wao -a PE_windows_fst01.txt -b "$gff3" | awk '{OFS="\t"; split($14, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$9,$10,$15,$12,$4,$5,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_fst_header.txt - > PE_intersect_fst01.txt
bedtools intersect -sorted -wao -a PE_windows_fst_all.txt -b "$gff3" | awk '{OFS="\t"; split($14, a, ";"); id=a[1]; split(id, b, "ID="); full=b[2]; split(full, c, "."); gene=c[1]; print $1,$2,$3,$9,$10,$15,$12,$4,$5,gene;}' - | sort -k1,1 -k2n,2 - | cat intersect_fst_header.txt - > PE_intersect_all.txt

# gene information is in the form: ID=ppa023343m.g;Name=ppa023343m.g;Dbxref="Phytozome:ppa023343m.g"
# http://www.unix.com/shell-programming-and-scripting/17696-split-field-awk-script.html

#intersect output with thetas file query		intersect output with fst file query
#col1	scaffold number					col1   scaffold number
#col2	winstar						col2   winstar
#col3	winend						col3   winend
#col4	thetaW						col4   unweighted
#col5	thetaPi						col5   weighted
#col6	Tajima						------------------------------------
#col7	fayh						col6   scaffold number
#col8	zeng                                            col7   genome
#--------------------------------------                 col8   gene
#col9	scaffold number                                 col9   start
#col10	genome                                          col10  stop
#col11	gene                                            col11  . (?)
#col12	start                                           col12  strand
#col13	stop						col13  .,0,1 (?)
#col14	. (?)						col14  gene ID, etc
#col15	strand						-----------------------------------
#col16	.,0,1 (?)					col15  intersect (bp)
#col17	gene ID, etc.
#--------------------------------------
#col18 intersect (bp)


#final file
#new	old	feature
#col1	col1   scaffold number
#col2	col2   winstar
#col3	col3   winend
#col4	col12  start
#col5	col13  stop
#col6	col18	intersect (bp)
#col7	col15  strand
#-----------------------------
#col8	col4   thetaW
#col9	col5   thetaPi
#col10	col6   Tajima
#col11	col7   fayh
#col12	col8   zeng
#-----------------------------
#col13	col17  gene ID, etc.

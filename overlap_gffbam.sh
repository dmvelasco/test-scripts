#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/genetree
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stdout-extractgenes.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%A_%a-stderr-extraxtgenes.txt
#SBATCH -J genext
#SBATCH -p bigmemm
#SBATCH -a 9
#SBATCH -t 1-00:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
set -e
set -u

# Declare number variables
x=$SLURM_ARRAY_TASK_ID
i=$(( x-1 ))


# Load modules
module load zlib

# Declare directories and file prefix
dir1="/home/dmvelasc/Software/bedtools2/bin"	# bedtools location
dir2="/home/dmvelasc/Projects/Prunus/Data/BAM"		# BAM directory
dir3="/home/dmvelasc/Data/references/persica-SCF"	# gff3 location

# Declare prefix array
#declare -a main=(DPRU0194 DPRU0579 DPRU0582 DPRU1467.9 DPRU1871.1 DPRU2327.16 DPRU2493.7 DPRU2578.2 Lovell fenzliana TNP arabica DPRU1791.3 DPRU2374.12 DPRU1456.4 DPRU2301 DPRU1462.2 DPRU1207.2 DPRU2331.9 DPRU0210)
declare -a abbr=(PR01 PC01 PS02 PK01 PU01 PT01 PV02 PD01 PP15 PF01 PD02 PB01 PD03 PD04 PD05 PD06 PD07 PD08 PD09 PD10)
# Public sequences for later
# Public sequences for later - SOME HAVE MULTIPLE SRA RUNS
#declare -a sra=(SRR765861 SRR765850 SRR765838 SRR765679 SRR502998 SRR502999 SRR502985 SRR502994 SRR502992 SRR502993 SRR502990 SRR502991 SRR502987 SRR502989 SRR502986 SRR503000 SRR503001 SRR502983 SRR502997 SRR502995 SRR502996 SRR501836 SRR068361 SRR068359 SRR068360 SRR502984 SRR502982)
#declare -a pub=(PD11 PD12 PD13 PD14 PG01.1 PG01.2 PP01 PP02 PP03.1 PP03.2 PP04.1 PP04.2 PP05.1 PP05.2 PP06 PP07.1 PP07.2 PP08 PP09 PP10.1 PP10.2 PP11 PP12 PP13 PP14 PS01 PV01)
# ALSO ADD PUBLIC P. MIRA, NEW P. DAVIDIANA SEQUENCES

acc="${abbr["$i"]}"

# general idea of how to extract genes from bam
# from http://bedtools.readthedocs.io/en/latest/content/overview.html#overlapping-intersecting-features
#basic command
#bedtools intersect -a alignedReads.bam -b exons.bed

#command for cds, once cds extracted
#bedtools intersect -a "$acc".bam -b exons.gff > out_exons.bed #use locally reassembled bam from GATK here

# command for whole gene intersection
# use locally reassembled bam from GATK HaplotypeCaller here
"$dir1"/bedtools intersect -a "$dir2"/"$acc"_sorted_markdup.bam -b "$dir3"/Prunus_persica_v1.0_genes_only.gff3 > "$acc"_out_genes.bam

# try interactively because farm queue impacted...
# below works, once finally figured out the copying changed the '-' in the command line
#/home/dmvelasc/Software/bedtools2/bin/bedtools intersect -a /home/dmvelasc/Projects/Prunus/Data/BAM/PP15_sorted_markdup.bam -b /home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_genes_only.gff3 > PP15_gene_overlap.bam

# if this works, then to fasta format...

# FASTA format
# samtools fasta PP15_gene_overlap.bam > PP15_gene_overlap.fa
# generates the following:
#>FCD2GPRACXX:7:1110:13468:92406#/1
#AAAAATTCCCTCATCGATTGCTTAAGTGGGTGTAGCAAAAAATTGCTTAACCGGGAGAGTCAGAGAGAGGGAGAGGTGTAAAGCAATGGC
#>FCD2GPRACXX:7:1203:18773:68439#/2
#TGGTAGGTCAGACGCGGAGAGCGCCATTGCTTTACACCTCTCCCTCTCTCTGACTCTCCCGGTTAAGCAATTTTTTGCTACACCCACTTA
#>FCD2GPRACXX:7:2204:17884:81149#/2
#AGCAAAAAATTGCTTAACCGGGAGAGTCAGAGAGAGGGAGAGGTGTAAAGCAATGGCGCTCTCCGCGTCTGACCTACCAGCCATGTACTC
#>FCD2GPRACXX:8:2106:3577:55953#/1
#GAGTACATGGCTGGTAGGTCAGACGCGGAGAGCGCCATTGCTTTACACCTCTCCCTCTCTCTGACTCTCCCGGTTAAGCAATTTTTTGCT
#>FCD2GPRACXX:7:1301:2533:13561#/1
#AGAGAGTACATGGCTGGTAGGTCAGACGCGGAGAGCGCCATTGCTTTACACCTCTCCCTCTCTCTGACTCTCCCGGTTAAGCAATTTT

# PERHAPS ANOTHER OPTION...
java -Xmx2g -jar GenomeAnalysisTK.jar \
-R MY_REFERENCE.fa \
-T FastaAlternateReferenceMaker \
-o MY_REFERENCE_WITH_SNPS_FROM_VCF.fa \
--variant MY_VCF_IN_VCF_4.0_FORMAT.vcf
# from: https://www.biostars.org/p/17705/
# still how to get specific regions..., -L option with regions from gff?
# does -L take region file? sort of: -T SelectVariants \ -L /path/to/my.interval_list
# from: http://gatkforums.broadinstitute.org/gatk/discussion/2441/can-selectvariants-be-used-to-limit-vcf-files-by-interval-list
# see: https://software.broadinstitute.org/gatk/documentation/article?id=1319

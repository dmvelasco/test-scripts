#!/bin/bash -l
#SBATCH -D /group/jrigrp3/Velasco/Prunus/fastq/
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-cat.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-cat.txt
#SBATCH -J orig
#SBATCH -p bigmemm
#SBATCH -t 3-00:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
set -e
set -u

mkdir -p /scratch/dmvelasc/
temp="/scratch/dmvelasc"

# concatenate original sequencing FASTQ files for each sample
# keeping forward (_1) and reverse (_2) separate
# move to a temporary directory (concat/) in working directory

#DPRU0194/
#gunzip -c DPRU0194/FCD2GPRACXX-SZAIPI032084-64_L3_1.fq.gz DPRU0194/FCD2GPRACXX-SZAIPI032084-64_L4_1.fq.gz | gzip > "$temp"/DPRU0194_1.fastq.gz
#mv "$temp"/DPRU0194_1.fastq.gz concat/
#gunzip -c DPRU0194/FCD2GPRACXX-SZAIPI032084-64_L3_2.fq.gz DPRU0194/FCD2GPRACXX-SZAIPI032084-64_L4_2.fq.gz | gzip > "$temp"/DPRU0194_2.fastq.gz
#mv "$temp"/DPRU0194_2.fastq.gz concat/

#DPRU0579/
#cp DPRU0579/FCC2GP0ACXX-SZAIPI032110-94_L8_1.fq.gz concat/DPRU0579_1.fastq.gz
#cp DPRU0579/FCC2GP0ACXX-SZAIPI032110-94_L8_2.fq.gz concat/DPRU0579_2.fastq.gz
# no need to concatenate, only one original sequence each for forward and reverse
# just copy as new file name

#DPRU0582/
#gunzip -c DPRU0582/FCD2H8PACXX-SZAIPI034950-14_L7_1.fq.gz DPRU0582/FCD2H8PACXX-SZAIPI034950-14_L8_1.fq.gz | gzip > "$temp"/DPRU0582_1.fastq.gz
#mv "$temp"/DPRU0582_1.fastq.gz concat/
#gunzip -c DPRU0582/FCD2H8PACXX-SZAIPI034950-14_L7_2.fq.gz DPRU0582/FCD2H8PACXX-SZAIPI034950-14_L8_2.fq.gz | gzip > "$temp"/DPRU0582_2.fastq.gz
#mv "$temp"/DPRU0582_2.fastq.gz concat/

#DPRU1467.9/ ######################## HAD A PROBLEM SAID NOT GZIP FORMAT - originally piped to gunzip instead of gzip :( ################################
#gunzip -c DPRU1467.9/FCD2GPRACXX-SZAIPI032087-79_L3_1.fq.gz DPRU1467.9/FCD2GPRACXX-SZAIPI032087-79_L4_1.fq.gz | gzip > "$temp"/DPRU1467.9_1.fastq.gz
#mv "$temp"/DPRU1467.9_1.fastq.gz concat/
#gunzip -c DPRU1467.9/FCD2GPRACXX-SZAIPI032087-79_L3_2.fq.gz DPRU1467.9/FCD2GPRACXX-SZAIPI032087-79_L4_2.fq.gz | gzip > "$temp"/DPRU1467.9_2.fastq.gz
#mv "$temp"/DPRU1467.9_2.fastq.gz concat/

#DPRU1871.1/
#gunzip -c DPRU1871.1/FCD2H8PACXX-SZAIPI034951-15_L7_1.fq.gz DPRU1871.1/FCD2H8PACXX-SZAIPI034951-15_L8_1.fq.gz | gzip > "$temp"/DPRU1871.1_1.fastq.gz
#mv "$temp"/DPRU1871.1_1.fastq.gz concat/
#gunzip -c DPRU1871.1/FCD2H8PACXX-SZAIPI034951-15_L7_2.fq.gz DPRU1871.1/FCD2H8PACXX-SZAIPI034951-15_L8_2.fq.gz | gzip > "$temp"/DPRU1871.1_2.fastq.gz
#mv "$temp"/DPRU1871.1_2.fastq.gz concat/

#DPRU2327.16/
#gunzip -c DPRU2327.16/FCD2GPRACXX-SZAIPI032085-66_L3_1.fq.gz DPRU2327.16/FCD2GPRACXX-SZAIPI032085-66_L4_1.fq.gz | gzip > "$temp"/DPRU2327.16_1.fastq.gz
#mv "$temp"/DPRU2327.16_1.fastq.gz concat/
#gunzip -c DPRU2327.16/FCD2GPRACXX-SZAIPI032085-66_L3_2.fq.gz DPRU2327.16/FCD2GPRACXX-SZAIPI032085-66_L4_2.fq.gz | gzip > "$temp"/DPRU2327.16_2.fastq.gz
#mv "$temp"/DPRU2327.16_2.fastq.gz concat/

#DPRU2493.7/
#gunzip -c DPRU2493.7/FCD2GPRACXX-SZAIPI032086-74_L3_1.fq.gz DPRU2493.7/FCD2GPRACXX-SZAIPI032086-74_L4_1.fq.gz | gzip > "$temp"/DPRU2493.7_1.fastq.gz
#mv "$temp"/DPRU2493.7_1.fastq.gz concat/
#gunzip -c DPRU2493.7/FCD2GPRACXX-SZAIPI032086-74_L3_2.fq.gz DPRU2493.7/FCD2GPRACXX-SZAIPI032086-74_L4_2.fq.gz | gzip > "$temp"/DPRU2493.7_2.fastq.gz
#mv "$temp"/DPRU2493.7_2.fastq.gz concat/

#DPRU2578.2/
#gunzip -c DPRU2578.2/FCD2H8PACXX-SZAIPI034949-13_L7_1.fq.gz DPRU2578.2/FCD2H8PACXX-SZAIPI034949-13_L8_1.fq.gz | gzip > "$temp"/DPRU2578.2_1.fastq.gz
#mv "$temp"/DPRU2578.2_1.fastq.gz concat/
#gunzip -c DPRU2578.2/FCD2H8PACXX-SZAIPI034949-13_L7_2.fq.gz DPRU2578.2/FCD2H8PACXX-SZAIPI034949-13_L8_2.fq.gz | gzip > "$temp"/DPRU2578.2_2.fastq.gz
#mv "$temp"/DPRU2578.2_2.fastq.gz concat/

#FPS-Lovell/
#gunzip -c FPS-Lovell/FCD2GPRACXX-SZAIPI032088-80_L7_1.fq.gz FPS-Lovell/FCD2GPRACXX-SZAIPI032088-80_L8_1.fq.gz | gzip > "$temp"/FPS-Lovell_1.fastq.gz
#mv "$temp"/FPS-Lovell_1.fastq.gz concat/
#gunzip -c FPS-Lovell/FCD2GPRACXX-SZAIPI032088-80_L7_2.fq.gz FPS-Lovell/FCD2GPRACXX-SZAIPI032088-80_L8_2.fq.gz | gzip > "$temp"/FPS-Lovell_2.fastq.gz
#mv "$temp"/FPS-Lovell_2.fastq.gz concat/

#Sample_JRIDV1_A
gunzip -c Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_001.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_002.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_003.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_004.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_005.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_006.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_007.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R1_008.fastq.gz | gzip > "$temp"/Sample_A_1.fastq.gz
mv "$temp"/Sample_A_1.fastq.gz concat/
gunzip -c Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_001.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_002.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_003.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_004.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_005.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_006.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_007.fastq.gz Sample_JRIDV1_A/JRIDV1_A_CGATGT_L005_R2_008.fastq.gz | gzip > "$temp"/Sample_A_2.fastq.gz
mv "$temp"/Sample_A_2.fastq.gz concat/

#Sample_JRIDV1_B
gunzip -c Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_001.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_002.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_003.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_004.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_005.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_006.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_007.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R1_008.fastq.gz | gzip > "$temp"/Sample_B_1.fastq.gz
mv "$temp"/Sample_B_1.fastq.gz concat/
gunzip -c Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_001.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_002.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_003.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_004.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_005.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_006.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_007.fastq.gz Sample_JRIDV1_B/JRIDV1_B_TGACCA_L005_R2_008.fastq.gz | gzip > "$temp"/Sample_B_2.fastq.gz
mv "$temp"/Sample_B_2.fastq.gz concat/

#Sample_JRIDV1_C/
#gunzip -c Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_001.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_002.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_003.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_004.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_005.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_006.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R1_007.fastq.gz | gzip > "$temp"/Sample_C_1.fastq.gz
#mv "$temp"/Sample_C_1.fastq.gz concat/
#gunzip -c Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_001.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_002.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_003.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_004.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_005.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_006.fastq.gz Sample_JRIDV1_C/JRIDV1_C_ACAGTG_L005_R2_007.fastq.gz | gzip > "$temp"/Sample_C_2.fastq.gz
#mv "$temp"/Sample_C_2.fastq.gz concat/

#Sample_JRIDV1_D/
#gunzip -c Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_001.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_002.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_003.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_004.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_005.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_006.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_007.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R1_008.fastq.gz | gzip > "$temp"/Sample_D_1.fastq.gz
#mv "$temp"/Sample_D_1.fastq.gz concat/
#gunzip -c Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_001.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_002.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_003.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_004.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_005.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_006.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_007.fastq.gz Sample_JRIDV1_D/JRIDV1_D_GCCAAT_L005_R2_008.fastq.gz | gzip > "$temp"/Sample_D_2.fastq.gz
#mv "$temp"/Sample_D_2.fastq.gz concat/

#Sample_JRIDV1_E/
#gunzip -c Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_001.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_002.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_003.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_004.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_005.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_006.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_007.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_008.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R1_009.fastq.gz | gzip > "$temp"/Sample_E_1.fastq.gz
#mv "$temp"/Sample_E_1.fastq.gz concat/
#gunzip -c Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_001.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_002.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_003.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_004.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_005.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_006.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_007.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_008.fastq.gz Sample_JRIDV1_E/JRIDV1_E_CAGATC_L005_R2_009.fastq.gz | gzip > "$temp"/Sample_E_2.fastq.gz
#mv "$temp"/Sample_E_2.fastq.gz concat/

#Sample_JRIDV1_F/
#gunzip -c Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_001.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_002.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_003.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_004.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_005.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_006.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R1_007.fastq.gz | gzip > "$temp"/Sample_F_1.fastq.gz
#mv "$temp"/Sample_F_1.fastq.gz concat/
#gunzip -c Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_001.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_002.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_003.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_004.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_005.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_006.fastq.gz Sample_JRIDV1_F/JRIDV1_F_CTTGTA_L005_R2_007.fastq.gz | gzip > "$temp"/Sample_F_2.fastq.gz
#mv "$temp"/Sample_F_2.fastq.gz concat/

#Sample_JRIDV1_G/
#gunzip -c Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_001.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_002.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_003.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_004.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_005.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_006.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_007.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R1_008.fastq.gz | gzip > "$temp"/Sample_G_1.fastq.gz
#mv "$temp"/Sample_G_1.fastq.gz concat/
#gunzip -c Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_001.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_002.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_003.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_004.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_005.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_006.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_007.fastq.gz Sample_JRIDV1_G/JRIDV1_G_AGTCAA_L005_R2_008.fastq.gz | gzip > "$temp"/Sample_G_2.fastq.gz
#mv "$temp"/Sample_G_2.fastq.gz concat/

#Sample_JRIDV1_H/
#gunzip -c Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_001.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_002.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_003.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_004.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_005.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_006.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_007.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_008.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_009.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R1_010.fastq.gz | gzip > "$temp"/Sample_H_1.fastq.gz
#mv "$temp"/Sample_H_1.fastq.gz concat/
#gunzip -c Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_001.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_002.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_003.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_004.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_005.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_006.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_007.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_008.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_009.fastq.gz Sample_JRIDV1_H/JRIDV1_H_AGTTCC_L005_R2_010.fastq.gz | gzip > "$temp"/Sample_H_2.fastq.gz
#mv "$temp"/Sample_H_2.fastq.gz concat/

#UCD-TNP/
#gunzip -c UCD-TNP/FCC2GP0ACXX-SZAIPI032108-96_L6_1.fq.gz UCD-TNP/FCC2GP0ACXX-SZAIPI032108-96_L7_1.fq.gz | gzip > "$temp"/TNP_1.fastq.gz
#mv "$temp"/TNP_1.fastq.gz concat/
#gunzip -c UCD-TNP/FCC2GP0ACXX-SZAIPI032108-96_L6_2.fq.gz UCD-TNP/FCC2GP0ACXX-SZAIPI032108-96_L7_2.fq.gz | gzip > "$temp"/TNP_2.fastq.gz
#mv "$temp"/TNP_2.fastq.gz concat/

#UCD-fenzliana/
#gunzip -c UCD-fenzliana/FCD2GPRACXX-SZAIPI032089-81_L3_1.fq.gz UCD-fenzliana/FCD2GPRACXX-SZAIPI032089-81_L4_1.fq.gz | gzip > "$temp"/fenzliana_1.fastq.gz
#mv "$temp"/fenzliana_1.fastq.gz concat/
#gunzip -c UCD-fenzliana/FCD2GPRACXX-SZAIPI032089-81_L3_2.fq.gz UCD-fenzliana/FCD2GPRACXX-SZAIPI032089-81_L4_2.fq.gz | gzip > "$temp"/fenzliana_2.fastq.gz
#mv "$temp"/fenzliana_2.fastq.gz concat/

#USDA-arabica/
#gunzip -c USDA-arabica/FCC2GP0ACXX-SZAIPI032109-95_L6_1.fq.gz USDA-arabica/FCC2GP0ACXX-SZAIPI032109-95_L7_1.fq.gz | gzip > "$temp"/arabica_1.fastq.gz
#mv "$temp"/arabica_1.fastq.gz concat/
#gunzip -c USDA-arabica/FCC2GP0ACXX-SZAIPI032109-95_L6_2.fq.gz USDA-arabica/FCC2GP0ACXX-SZAIPI032109-95_L7_2.fq.gz | gzip > "$temp"/arabica_2.fastq.gz
#mv "$temp"/arabica_2.fastq.gz concat/

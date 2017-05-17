#!/bin/bash -l
#SBATCH -D /group/jrigrp7/josh/data/bams
#SBATCH -o /group/jrigrp7/josh/logs/teo_chr_1_stdout_%A_%a.txt
#SBATCH -e /group/jrigrp7/josh/logs/teo_chr_1_stderr_%A_%a.txt
#SBATCH -J chr_1
#SBATCH --array=1-5
#SBATCH --cpus-per-task=3

set -e
set -u
set -o pipefail

# set directories for reference genome, bam files, output, and scripts
# ref_dir=/home/jhough/projects/peach/data/reference
bam_dir=/group/jrigrp7/josh/data/bams
out_dir=/group/jrigrp7/josh/results
ref=/group/jrigrp7/josh/data/reference
# bamCaller=/home/jhough/projects/peach/msmc-tools/./bamCaller.py

################################################################################
fastq=$(find $fastq_dir -name "*.fastq" | xargs -n1 -I{} basename {} .fastq.gz | sort -s | sed -n "$SLURM_ARRAY_TASK_ID"p )

samtools mpileup -q 20 -Q 20 -C 50 -u -r scaffold_1 -f $ref ${fastq}.bam | \
bcftools call -c -V indels | vcffilter -f "DP > 15" | ${bamCaller} ${chr_1_depth} ${out_dir}/${peach_bam}_chr_1_mask.bed.gz | gzip -c > ${out_dir}/${peach_bam}_chr_1.vcf.gz
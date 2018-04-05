#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemh
#SBATCH -t 4:00:00
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=64G
set -e
set -u

########## WRITTEN BY D. VELASCO ###########
## initial set up borrowed from M Stetter ##

####################
### Load modules ###
####################
# Anaconda 3, automatically enables SMC++ (and Mafft)
module load conda3
module load vcftools
module load tabix

#############################
### Set up the parameters ###
#############################

####### ARRAYS #######

##############################
##### ACQUIRE SETUP DATA #####
## TWO FOR SPLIT ##
# path to sample list
list="/home/dmvelasc/Projects/Prunus/Script/smcpp_data.txt"
# view list to select pops/subpops, important below
# PD	sub1-10			lines 1-10
# PD	sub1	range		line 1
# PD	sub7	China		line 7
# PD	sub8	C Asia		line 8
# PD	sub10	S Europe/W Asia	line 10
# PP	sub1-7			lines 11-17
# PP	sub1	range		line 11
# PP	sub3	China		line 13
# PP	sub4	China/Korea	line 14
# PM	sub1-3			lines 18-20
# PV	sub1-2			lines 21-22
# PS	all			line 23
# PG	all			line 24

line1=1		#line number of first pop/subpop wanted
line2=21	#line number of second pop/subpop wanted

# mapfile to extract first pop/subpop information and send to an array
mapfile -s "$line1" -n 1 -t id < "${list}"
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop1="${arr[0]}"
sub1="${arr[1]}"
sample_1a="${arr[2]}"
sample_1b="${arr[3]}"
sample_1c="${arr[4]}"
sample_1d="${arr[5]}"

# mapfile to extract second pop/subpop information and send to an array
mapfile -s "$line2" -n 1 -t id < "${list}"
arr=(`echo "${id[0]}"`)

# declare variables, created from array
pop2="${arr[0]}"
sub2="${arr[1]}"
sample_2a="${arr[2]}"
sample_2b="${arr[3]}"
sample_2c="${arr[4]}"
sample_2d="${arr[5]}"
################################

####### PATHS #######
# genome reference file
genome="/home/dmvelasc/Data/references/persica-SCF/Prunus_persica_v1.0_scaffolds.fa"
# Joint VCF file - test joint VCF file
vcf="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_jointcalls.vcf"
# filtered joint VCF file - dulcis test VCF file
vcf_filt="/home/dmvelasc/Projects/Prunus/Analysis/VCF_GATK/test_dulcis.recode.vcf"
# SMC++ prepped file
smc_in="/home/dmvelasc/Projects/Prunus/Data/smcpp_input/"
final="mu1_38x-8"

####### PARAMETERS #######
mu="7.77e-9"		# population mutation rate
cut="5000"		# cutoff length for homozygosity

####################
### Begin script ###
####################
echo -e "begin SMC++ preparation\n get individuals"
date

# VCF filter for pairwise comparison
sub="all" #subset name
vcftools --vcf "$vcf" \
--indv "$sample_1a" \
--indv "$sample_1b" \
--indv "$sample_1c" \
--indv "$sample_1d" \
--indv "$sample_2a" \
--indv "$sample_2b" \
--indv "$sample_2c" \
--indv "$sample_2d" \
--min-alleles 2 --max-alleles 2 \
--recode \
--out split-"$pop1"_"$sub1"-"$pop2"_"$sub2"

##### NEEDED FOR INITIAL PREP #####
mv /home/dmvelasc/Projects/Prunus/Analysis/smcpp/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf "$vcf_filt"/


echo -e "convert vcf file to SMC++ format file"
date

bgzip -f "$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf > "$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf.gz
tabix -fp vcf "$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf.gz

echo -e "Run for loop by chromosome as per smcpp instructions"
date

##### NEEDED FOR INITIAL PREP #####
# select correct configuration according to paired populations
for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff "$cut" \
"$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf.gz \
"$smc_in"/split-"$pop1"_"$sub1"-"$i".smc.gz scaffold_"$i" \
"$pop1":"$sample_1a","$sample_1b","$sample_1c","$sample_1d"
  smc++ vcf2smc --missing-cutoff "$cut" \
"$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf.gz \
"$smc_in"/split-"$pop2"_"$sub2"-"$i".smc.gz scaffold_"$i" \
"$pop2":"$sample_2a","$sample_2b","$sample_2c","$sample_2d"

###########
# SEE BELOW
#  smc++ estimate -o "$pop1"_"$sub1" "$mu" "$smc_in"/split-"$pop1"_"$sub1"-"$i".smc.gz
#  smc++ estimate -o "$pop2"_"$sub2" "$mu" "$smc_in"/split-"$pop2"_"$sub2"-"$i".smc.gz
###########
done


echo -e "begin SMC++ analysis"
date
# SMC++ analysis
smc++ estimate -o smc_analysis/"$pop1"_"$sub1" "$mu" "$smc_in"/split-"$pop1"_"$sub1"-*.smc.gz
smc++ estimate -o smc_analysis/"$pop2"_"$sub2" "$mu" "$smc_in"/split-"$pop2"_"$sub2"-*.smc.gz

#--polarization-error 0.5
# --polarization-error: if the identity of the ancestral allele is not known,
# these options can be used to specify a prior over it. With polarization error p,
# emissions probabilities for entry CSFS(a,b) will be computed as
# (1-p) CSFS(a,b) + p CSFS(2-a, n-b). The default setting is 0.5,
# i.e. the identity of the ancestral allele is not known.
# --unfold is an alias for --polarization-error 0. If the ancestral allele is known
# (from an outgroup, say) then this option will use the unfolded SFS for computing
# probabilities. Incorrect usage of this feature may lead to erroneous results.
# $mu is per generation mutation rate, will probably need to run with three different values based on Xie et al.

################## SPLIT #####################
# substitute different individuals for different combinations
for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff "$cut" \
"$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf.gz \
"$smc_in"/split_"$pop1"_"$sub1"-"$pop2"_"$sub2"-"$i".smc.gz scaffold_"$i" \
"$pop1":"$sample_1a","$sample_1b","$sample_1c","$sample_1d" \
"$pop2":"$sample_2a","$sample_2b","$sample_2c","$sample_2d"
  smc++ vcf2smc --missing-cutoff "$cut" \
"$vcf_filt"/split-"$pop1"_"$sub1"-"$pop2"_"$sub2".recode.vcf.gz \
"$smc_in"/split_"$pop2"_"$sub2"-"$pop1"_"$sub1"-"$i".smc.gz scaffold_"$i" \
"$pop2":"$sample_2a","$sample_2b","$sample_2c","$sample_2d" \
"$pop1":"$sample_1a","$sample_1b","$sample_1c","$sample_1d"
done

# Run split to refine the marginal estimates
smc++ split -o split/ smc_analysis/"$pop1"_"$sub1"/model.final.json smc_analysis/"$pop2"_"$sub2"/model.final.json "$smc_in"/*.smc.gz
smc++ plot -c joint.pdf split/model.final.json

##### posterior #####
# export (and visualize) the posterior distributon of the TMRCA

#smc++ posterior [--unfold | --polarization-error p] [--start START] [--end END]
# [--thinning k] [--heatmap heatmap.(pdf|png|gif|jpeg)] [--colorbar] [--M M]
# model.final.json arrays.npz data.smc[.gz] [data.smc[.gz] ...]

#positional arguments:
#  model.final.json      SMC++ model to use in forward-backward algorithm
#  arrays.npz            location to save posterior decoding arrays
#  data.smc[.gz]         SMC++ data set(s) to decode

#optional arguments:
#  -h, --help            show this help message and exit
#  -v, --verbose         increase debugging output, specify multiply times for more
#  --unfold              use unfolded SFS (alias for -p 0.0)
#  --polarization-error p, -p p
#                        uncertainty parameter for polarized SFS: observation (a,b) has probability [(1-p)*CSFS_{a,b} +
#                        p*CSFS_{2-a,n-2-b}]. default: 0.5
#  --start START         base at which to begin posterior decode
#  --end END             base at which to end posterior decode
#  --thinning k          emit full SFS only every <k>th site. default: 1
#  --heatmap heatmap.(pdf|png|gif|jpeg)
#                        Also draw a heatmap of the posterior TMRCA.
#  --colorbar            If plotting, add a colorbar

#HMM parameters:
#  --M M                 number of hidden states

################## SPLIT #####################


echo -e "plot SMC++ results"
date

# move files output files to subdirectory
mkdir -p smc_analysis/"$final"
mv smc_analysis/model.final.json smc_analysis/"$final"
mv smc_analysis/.model.iter*.json smc_analysis/"$final"
mv smc_analysis/.debug.txt smc_analysis/"$final"

#-g	sets generation time in years to scale x-axis, otherwise in coalescent units
#--logy	plots the y-axis on a log scale
#-c	produces CSV-formatted table containing the data used to generate the plot


################################################
################################################

### M STETTER SCRIPT BELOW


#!/bin/bash -l
#SBATCH -D /home/mstetter/amaranth_domestication/
#SBATCH -o /home/mstetter/amaranth_domestication/data/smcpp/stdout-%j.txt
#SBATCH -J run_smcpp
#SBATCH -t 120:00:00
#SBATCH --mem 50gb
set -e
set -u

#module load tabix
# run below command on farm before starting the script
#source /home/mstetter/tools/smcpp/bin/activate
# file must be  bgziped and indexted with tabix

#gunzip data/wgrs_snps/16chr_ann_gatk_filter_biallelicSNPs_noINDELs.vcf.gz
#bgzip data/wgrs_snps/16chr_ann_gatk_filter_biallelicSNPs_noINDELs.vcf
#tabix -fp vcf data/wgrs_snps/16chr_ann_gatk_filter_biallelicSNPs_noINDELs.vcf.gz

chr_length=(38124660 35657244 30204323 28349311 25672467 24628041 24364990 23766980 22691259 22670516 22280117 22052327 20679869 20190685 17522127 16951160)
populations=(caudatus cruentus hypochondriacus quitensis hybridus)

OUTPUTfolder=data/smcpp/
VCFfile=data/wgrs_snps/16chr_ann_gatk_filter_biallelicSNPs_noINDELs.vcf.gz
MISSING=1000
# DO this part only once if needed
#mkdir -p data/population_subsamples
#population_files=(caudatus_samples cruentus_samples_corrected hypochondriacus_samples_corrected quitensis_samples hybridus_samples)

#for i in {0..4};do
#shuf -n 19 data/population_files/${population_files["$i"]} > data/population_subsamples/${populations["$i"]}
#done

mkdir -p $OUTPUTfolder/figures_output
mkdir -p $OUTPUTfolder/smc_input
mkdir -p $OUTPUTfolder/smc_output

for i in "${populations[@]}";do
  SAMPLES=($(cat data/population_subsamples/"$i"))
  SAMPLES=$(IFS=, ; echo "${SAMPLES[*]}")
  echo $SAMPLES

  for CHR in {1..16};do
    smc++ vcf2smc --missing-cutoff $MISSING --length ${chr_length[$CHR-1]} $VCFfile $OUTPUTfolder/smc_input/"$i"_chr"$CHR".smc.gz $CHR "$i":"$SAMPLES"
    echo Chromosome $CHR finished
  done

  mkdir -p $OUTPUTfolder/smc_output/"$i"/
  smc++ estimate -o $OUTPUTfolder/smc_output/"$i"/ 7e-9 $OUTPUTfolder/smc_input/"$i"_*.smc.gz
  smc++ plot -c $OUTPUTfolder/figures_output/"$i"_plot.png $OUTPUTfolder/smc_output/"$i"/model*.json
done

##### Calculate splits

HYBR=($(cat data/population_subsamples/hybridus))
HYBR=$(IFS=, ; echo "${HYBR[*]}")

mkdir -p $OUTPUTfolder/smc_input/splits/

for SPLITPOP in "${populations[@]:0:4}";do
  echo $SPLITPOP

  SAMPLES=($(cat data/population_subsamples/"$SPLITPOP"))
  SAMPLES=$(IFS=, ; echo "${SAMPLES[*]}")
  echo $SAMPLES
  mkdir -p $OUTPUTfolder/smc_input/splits/split_hybridus_"$SPLITPOP"

  for CHR in {1..16};do
    echo $CHR
    smc++ vcf2smc --missing-cutoff $MISSING --length ${chr_length[$CHR-1]} $VCFfile $OUTPUTfolder/smc_input/splits/split_hybridus_"$SPLITPOP"/hybridus_"$SPLITPOP"_chr"$CHR".smc.gz "$CHR" hybridus:"$HYBR" "$SPLITPOP":"$SAMPLES"
    smc++ vcf2smc --missing-cutoff $MISSING --length ${chr_length[$CHR-1]} $VCFfile $OUTPUTfolder/smc_input/splits/split_hybridus_"$SPLITPOP"/"$SPLITPOP"_hybridus_chr"$CHR".smc.gz "$CHR" "$SPLITPOP":"$SAMPLES" hybridus:"$HYBR"
  done

  mkdir -p $OUTPUTfolder/smc_output/split_"$SPLITPOP"_hyb/
  smc++ split -o $OUTPUTfolder/smc_output/split_"$SPLITPOP"_hyb/ $OUTPUTfolder/smc_output/hybridus/model.final.json $OUTPUTfolder/smc_output/"$SPLITPOP"/model.final.json $OUTPUTfolder/smc_input/splits/split_hybridus_"$SPLITPOP"/*.smc.gz
  smc++ plot -c $OUTPUTfolder/figures_output/split_"$SPLITPOP"_hybridus.pdf $OUTPUTfolder/smc_output/split_"$SPLITPOP"_hyb/model.final.json
done

## Plot everything
smc++ plot -c $OUTPUTfolder/figures_output/allinone.png $OUTPUTfolder/smc_output/*/model*.json

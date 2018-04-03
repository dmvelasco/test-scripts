#!/bin/bash -l
#SBATCH -D /home/dmvelasc/Projects/Prunus/Analysis/smcpp
#SBATCH -o /home/dmvelasc/Projects/Prunus/slurm-log/%j-stdout-smcpp.txt
#SBATCH -e /home/dmvelasc/Projects/Prunus/slurm-log/%j-stderr-smcpp.txt
#SBATCH -J smcpp
#SBATCH -p bigmemh
#SBATCH -t 4:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mail-user=dmvelasco@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=56G
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
# population array
declare -a pop=(PD PP PM PV PS PG)

# samples for each population
# PD, pop[0]
dulcis="PD02,PD03,PD04,PD05,PD06,PD07,PD08,PD09,PD10,PD11,PD12,PD13,PD14,PD16,PD17,PD18,PD20,PD21"
# PP, pop[1]
persica="PP02,PP03,PP04,PP05,PP06,PP08,PP11,PP13,PP14,PP15,PP37,PP38,PP39,PP40"
# PM, pop[2]
mira="PM01,PM02,PM03,PM04,PM05,PM06"
# PV, pop[3]
davidiana="PV01,PV02,PV03,PV04,PV05,PV06"
# PS, pop[4]
kansuensis="PS01,PS02,PS03,PS04"
# PG, pop[5]
ferganensis="PG02,PG03,PG04,PG05"

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
mu="7.77e-9"	# population mutation rate
cut="5000"	# cutoff length for homozygosity
pop1="${pop[2]}"	# population 1, P. mira
pop2="${pop[4]}"	# population 2, P. kansuensis

####################
### Begin script ###
####################
echo -e "begin SMC++ preparation\n get individuals"
date

# select individuals in pairwise combinations

# PD-PP; pairwise dulcis-persica; subset=all; 32 individuals; 30 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PD02 --indv PD03 --indv PD04 --indv PD05 --indv PD06 \
#--indv PD07 --indv PD08 --indv PD09 --indv PD10 --indv PD11 --indv PD12 --indv PD13 \
#--indv PD14 --indv PD16 --indv PD17 --indv PD18 --indv PD20 --indv PD21 --indv PP02 \
#--indv PP03 --indv PP04 --indv PP05 --indv PP06 --indv PP08 --indv PP11 --indv PP13 \
#--indv PP14 --indv PP15 --indv PP37 --indv PP38 --indv PP39 --indv PP40 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[0]}"-"${pop[1]}"

# PD-PV; pairwise dulcis-davidiana; subset=all; 24 individuals; 22 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PD02 --indv PD03 --indv PD04 --indv PD05 --indv PD06 \
#--indv PD07 --indv PD08 --indv PD09 --indv PD10 --indv PD11 --indv PD12 --indv PD13 \
#--indv PD14 --indv PD16 --indv PD17 --indv PD18 --indv PD20 --indv PV01 --indv PV02 \
#--indv PV03 --indv PV04 --indv PV05 --indv PV06 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[0]}"-"${pop[3]}"

# PP-PM; pairwise persica-mira; subset=all; 20 individuals; 18 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PP02 --indv PP03 --indv PP04 --indv PP05 --indv PP06 \
#--indv PP08 --indv PP11 --indv PP13 --indv PP14 --indv PP15 --indv PP37 --indv PP38 \
#--indv PP39 --indv PP40 --indv PM01 --indv PM02 --indv PM03 --indv PM04 --indv PM05 \
#--indv PM06 --min-alleles 2 --max-alleles 2 --recode --out "$sub"_"{pop[1]}"-"${pop[2]}"

# PP-PV; pairwise persica-davidiana; subset=all; 20 individuals; 18 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PP02 --indv PP03 --indv PP04 --indv PP05 --indv PP06 \
#--indv PP08 --indv PP11 --indv PP13 --indv PP14 --indv PP15 --indv PP37 --indv PP38 \
#--indv PP39 --indv PP40 --indv PV01 --indv PV02 --indv PV03 --indv PV04 --indv PV05 \
#--indv PV06 --min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[1]}"-"${pop[3]}"

# PP-PS; pairwise persica-kansuensis; subset=all; 18 individuals; 18 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PP02 --indv PP03 --indv PP04 --indv PP05 --indv PP06 \
#--indv PP08 --indv PP11 --indv PP13 --indv PP14 --indv PP15 --indv PP37 --indv PP38 \
#--indv PP39 --indv PP40 --indv PS01 --indv PS02 --indv PS03 --indv PS04 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"{pop[1]}"-"${pop[4]}"

# PP-PG; pairwise persica-ferganensis; subset=all; 18 individuals; 18 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PP02 --indv PP03 --indv PP04 --indv PP05 --indv PP06 \
#--indv PP08 --indv PP11 --indv PP13 --indv PP14 --indv PP15 --indv PP37 --indv PP38 \
#--indv PP39 --indv PP40 --indv PG02 --indv PG03 --indv PG04 --indv PG05 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[1]}"-"${pop[5]}"

# PM-PV; pairwise mira-davidiana; subset=all; 12 individuals; 12 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PM01 --indv PM02 --indv PM03 --indv PM04 --indv PM05 \
#--indv PM06 --indv PV01 --indv PV02 --indv PV03 --indv PV04 --indv PV05 --indv PV06 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[2]}"-"${pop[3]}"

# PM-PS; pairwise miras-kansuensis; subset=all; 10  individuals; 8 CPU?
sub="all" #subset name
vcftools --vcf "$vcf" --indv PM01 --indv PM02 --indv PM03 --indv PM04 --indv PM05 \
--indv PM06 --indv PS01 --indv PS02 --indv PS03 --indv PS04 \
--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[2]}"-"${pop[4]}"

# PM-PG; pairwise mira-ferganensis; subset=all; 10 individuals; 8 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PM01 --indv PM02 --indv PM03 --indv PM04 --indv PM05 \
#--indv PM06 --indv PG02 --indv PG03 --indv PG04 --indv PG05 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[2]}"-"${pop[5]}"

# PS-PV; pairwise kansuensis-davidiana; subset=all; 10  individuals; 8 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PS01 --indv PS02 --indv PS03 --indv PS04 --indv PV01 \
#--indv PV02 --indv PV03 --indv PV04 --indv PV05 --indv PV06 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[4]}"-"${pop[3]}"

# PS-PG; pairwise kansuensis-ferganensis; subset=all; 8 individuals; 6 CPU?
#sub="all" #subset name
#vcftools --vcf "$vcf" --indv PS01 --indv PS02 --indv PS03 --indv PS04 --indv PG02 \
#--indv PG03 --indv PG04 --indv PG05 \
#--min-alleles 2 --max-alleles 2 --recode --out "$sub"_"${pop[2]}"-"${pop[5]}"

##### NEEDED FOR INITIAL PREP #####
mv /home/dmvelasc/Projects/Prunus/Analysis/smcpp/"$sub"_"$pop1"-"$pop2".recode.vcf "$vcf_filt"/


echo -e "convert vcf file to SMC++ format file"
date

bgzip -f "$vcf_filt"/all_"$pop1"-"$pop2".recode.vcf > "$vcf_filt"/all_"$pop1".recode.vcf.gz
tabix -fp vcf "$vcf_filt"/all_"$pop1"-"$pop2".recode.vcf.gz

echo -e "Run for loop by chromosome as per smcpp instructions"
date

##### NEEDED FOR INITIAL PREP #####
# select correct configuration according to paired populations
for i in {1..8}; do
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/all_"$pop1"-"$pop2".recode.vcf.gz "$smc_in"/all_"$pop1"_"$i".smc.gz scaffold_"$i" "$pop1":"$mira"
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/all_"$pop1"-"$pop2".recode.vcf.gz "$smc_in"/all_"$pop2"_"$i".smc.gz scaffold_"$i" "pop2":"$kansuensis"
done


echo -e "begin SMC++ analysis"
date
# SMC++ analysis
smc++ estimate -o smc_analysis/"$pop1" "$mu" "$smc_in"/all_"$pop1"*.smc.gz
smc++ estimate -o smc_analysis/"$pop2" "$mu" "$smc_in"/all_"$pop2"*.smc.gz

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
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/all_"$pop1"-"$pop2".recode.vcf.gz \
"$smc_in"/all_"$pop1"-"$pop2"_"$i".smc.gz scaffold_"$i" "$pop1":"$mira" "$pop2":"$kansuensis"
  smc++ vcf2smc --missing-cutoff "$cut" "$vcf_filt"/all_"$pop1"-"$pop2".recode.vcf.gz \
"$smc_in"/all_"$pop2"-"$pop1"_"$i".smc.gz scaffold_"$i" "$pop2":"$kansuensis" "$pop1":"$mira"
done

# Run split to refine the marginal estimates
smc++ split -o "$pop1"-"$pop2"_split/ "$pop1"/model.final.json "$pop2"/model.final.json data/*.smc.gz
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

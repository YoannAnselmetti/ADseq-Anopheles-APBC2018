#!/bin/bash
###
###   Goal:
###      Script to run all scripts for Validation of ADseq algo with BLASTn align
###
###   INPUT:
###      1- Species name
###         (Anopheles_albimanus / Anopheles_arabiensis / Anopheles_dirus)
###      2- Species TAG
###         (Aalb / Aara / Adir)
###      3- Reads sampling TAG_ID
###         (50pourc / ALL)
###      4- Output directory where stats file will be stored 
###         (results/validation_ADseq/spi_20/stats/minia/with_ori)
###      5- AGP file 
###         (data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/50pourc/Anopheles_albimanus/BESST_output/pass2/info-pass2.agp)
###      6-  spi directory (Scaffolding Propagation Index)
###         (spi_default / spi_20)
###      7-  Boolean if minia contigs are original (minia) or merged (blastn)
###         (minia / blastn)
###      8-  Boolean if take into account orientation or not
###         (with_ori / without_ori)
###
###   OUTPUT:
###      - BLASTn results of CTG from MINIA assembly mapped on CTG from REF assembly to the gene ID
###
###   Name: run_for_1exp.sh                          Author: Yoann Anselmetti
###   Creation date: 2016/11/17                      Last modification: 2017/10/25
###



SCRIPTdir="bin/scripts/validation_ADseq/stats/01-compute_stats"
ADSEQdir="data/validation_ADseq/ADseq/BLASTn"
annotGfileORI=$ADSEQdir"/"$1"/"$3"/annotG_"$2"_ORI_filt"

echo $annotGfileORI

resultsADseqLinearized01="results/validation_ADseq/"$6"/"$2"/"$3"/ADseq/"$7"/validation_ADseq_"$2"_"$3"_Xtopo_Boltz_0.1_01_M1_kept"
resultsADseqLinearized05="results/validation_ADseq/"$6"/"$2"/"$3"/ADseq/"$7"/validation_ADseq_"$2"_"$3"_Xtopo_Boltz_0.1_05_M1_kept"
resultsADseqLinearized08="results/validation_ADseq/"$6"/"$2"/"$3"/ADseq/"$7"/validation_ADseq_"$2"_"$3"_Xtopo_Boltz_0.1_08_M1_kept"
# resultsADseqLinearized095="results/validation_ADseq/"$6"/"$2"/"$3"/ADseq/validation_ADseq_"$2"_"$3"_Xtopo_Boltz_0.1_M1_kept"
resultsADseqRAW="results/validation_ADseq/"$6"/"$2"/"$3"/ADseq/"$7"/validation_ADseq_"$2"_"$3"_Xtopo_Boltz_0.1.adjacencies.txt"

resultsADLinearized01="results/validation_ADseq/"$6"/"$2"/"$3"/ARt-DeCo/validation_ADseq_"$2"_"$3"_Xtopo_ARt-DeCo_Boltz_0.1_01_M1_kept"
resultsADLinearized05="results/validation_ADseq/"$6"/"$2"/"$3"/ARt-DeCo/validation_ADseq_"$2"_"$3"_Xtopo_ARt-DeCo_Boltz_0.1_05_M1_kept"
resultsADLinearized08="results/validation_ADseq/"$6"/"$2"/"$3"/ARt-DeCo/validation_ADseq_"$2"_"$3"_Xtopo_ARt-DeCo_Boltz_0.1_08_M1_kept"
# resultsADLinearized095="results/validation_ADseq/"$6"/"$2"/"$3"/ARt-DeCo/validation_ADseq_"$2"_"$3"_Xtopo_ARt-DeCo_Boltz_0.1_M1_kept"
resultsADRAW="results/validation_ADseq/"$6"/"$2"/"$3"/ARt-DeCo/validation_ADseq_"$2"_"$3"_Xtopo_ARt-DeCo_Boltz_0.1.adjacencies.txt"

resultsBESST="data/validation_ADseq/ADseq/BLASTn/"$1"/"$3"/scaff_"$2"_"$3"_3_"$7
outputAD=$4"/"$2"_"$3"_AD"
outputADseq=$4"/"$2"_"$3"_ADseq"
outputBESST=$4"/"$2"_"$3"_BESST"

mkdir -p $4

CTGfile="data/validation_ADseq/ADseq/BLASTn/"$1"/"$3"/CTG_"$2"_"$3"_filt"
ADJfileBESST=""
assoCTGfile="data/validation_ADseq/ADseq/BLASTn/"$1"/"$3"/assoCTG_"$2"_"$3

header="filter\t#PRED\t#FIND\t#FN\t#TP\t#FP\t#CFP\tCFP/FP\tRec\tPrec\tRecCFP\tPrecCFP"
echo -e $header > $outputAD
echo -e $header > $outputADseq
echo -e $header > $outputBESST

boolOri=""
if [ "$8" = "with_ori" ]; then
	boolOri="t"
elif [ "$8" = "without_ori" ]; then
	boolOri="f"
else
	echo "Parameter 8 should be equal to \"with_ori\" or \"without_ori\" and NOT: \""$8"\" !!!"
	exit
fi


#####
### MAIN
#####

# ADseq
$SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADseqLinearized01 $1 $4 0.1 $2 $3 ADseq $boolOri
$SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADseqLinearized05 $1 $4 0.5 $2 $3 ADseq $boolOri
$SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADseqLinearized08 $1 $4 0.8 $2 $3 ADseq $boolOri
# $SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADseqLinearized095 $1 $4 0.95 $2 $3 ADseq $boolOri
$SCRIPTdir/compare_newADJminia_ADJori_without_linearization.py $annotGfileORI $resultsADseqRAW $1 $outputADseq $boolOri

# ARt-DeCo
$SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADLinearized01 $1 $4 0.1 $2 $3 AD $boolOri
$SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADLinearized05 $1 $4 0.5 $2 $3 AD $boolOri
$SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADLinearized08 $1 $4 0.8 $2 $3 AD $boolOri
# $SCRIPTdir/compare_newADJminia_ADJori_with_linearization.py $annotGfileORI $resultsADLinearized095 $1 $$4 0.95 $2 $3 AD $boolOri
$SCRIPTdir/compare_newADJminia_ADJori_without_linearization.py $annotGfileORI $resultsADRAW $1 $outputAD $boolOri

# BESST
if [ "$7" = "minia" ]; then
	ADJfileBESST="data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/"$3"/"$1"/newADJ_BESST_complete"
	echo "$SCRIPTdir/transform_AGPfile_in_ADJfile_minia.py $5 $CTGfile $ADJfileBESST $assoCTGfile $1"
	$SCRIPTdir/transform_AGPfile_in_ADJfile_minia.py $5 $CTGfile $ADJfileBESST $assoCTGfile $1
elif [ "$7" = "blastn" ]; then
	ADJfileBESST="data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_k50/TRIMMOMATIC3/blastn/"$3"/"$1"/newADJ_BESST_complete"
	echo "$SCRIPTdir/transform_AGPfile_in_ADJfile_blastn.py $5 $CTGfile $ADJfileBESST $1"
	$SCRIPTdir/transform_AGPfile_in_ADJfile_blastn.py $5 $CTGfile $ADJfileBESST $1
else
	echo "ERROR, parameter 7 should be equal to \"minia\" or \"blastn\"!!!"
	exit
fi
echo "$SCRIPTdir/compare_newADJbesst_ADJori.py $annotGfileORI $resultsADseqLinearized01 $1 $ADJfileBESST $resultsBESST $2 $3 $4 $boolOri"
$SCRIPTdir/compare_newADJbesst_ADJori.py $annotGfileORI $resultsADseqLinearized01 $1 $ADJfileBESST $resultsBESST $2 $3 $4 $boolOri



#!/bin/bash
###
###   Goal:
###      Script to run all plot of Precision & Recall for the validation of ADseq algo
###
###   OUTPUT:
###      - ALL Precision & Recall plots
###
###   Name: run_all_exp.sh                          Author: Yoann Anselmetti
###   Creation date: 2017/01/27                     Last modification: 2017/10/25
###



SCRIPTdir="bin/scripts/validation_ADseq/stats/02-plot_precision_recall"
script1=$SCRIPTdir"/plot_Precision_Recall.py"
script2=$SCRIPTdir"/plot_Precision_ON_Recall.py"




boolORI="with_ori without_ori"
assembly="blastn"

for ori in $boolORI; do
	for ass in $assembly; do

		inputDIR="results/validation_ADseq/spi_20/stats/"$ass"/"$ori
		outputDIR="figures/precision_recall/spi_20/"$ass"/"$ori

		mkdir -p $inputDIR
		mkdir -p $outputDIR

		# ###################
		# ### PLOT Recall & Precision stats distribution
		# ###################
		echo "PLOT distribution of Recall & Precision stats:"
		echo -e "\tFor species Aalb:"
		echo -e "\t\tFor Scaffolding Propagation Index = 20:"
		echo -e "\t\t\t$script1 Anopheles_albimanus Aalb 50pourc $inputDIR $outputDIR"
		$script1 Anopheles_albimanus Aalb 50pourc $inputDIR $outputDIR

		echo -e "\t\t\t$script1 Anopheles_albimanus Aalb ALL $inputDIR $outputDIR"
		$script1 Anopheles_albimanus Aalb ALL $inputDIR $outputDIR



		echo -e "\n\tFor species Aara:"
		echo -e "\t\tFor Scaffolding Propagation Index = 20:"
		echo -e "\t\t\t$script1 Anopheles_arabiensis Aara 50pourc $inputDIR $outputDIR"
		$script1 Anopheles_arabiensis Aara 50pourc $inputDIR $outputDIR

		echo -e "\t\t\t$script1 Anopheles_arabiensis Aara ALL $inputDIR $outputDIR"
		$script1 Anopheles_arabiensis Aara ALL $inputDIR $outputDIR



		echo -e "\n\tFor species Adir:"
		echo -e "\t\tFor Scaffolding Propagation Index = 20:"
		echo -e "\t\t\t$script1 Anopheles_dirus Adir 50pourc $inputDIR $outputDIR:"
		$script1 Anopheles_dirus Adir 50pourc $inputDIR $outputDIR

		echo -e "\t\t\t$script1 Anopheles_dirus Adir ALL $inputDIR $outputDIR:"
		$script1 Anopheles_dirus Adir ALL $inputDIR $outputDIR



		###################
		### PLOT Precision/Recall
		###################
		echo ""
		echo "PLOT Precision/Recall stats:"
		echo -e "\tFor Scaffolding Propagation Index = 20:"
		echo -e "\t\t$script2 50pourc $inputDIR $outputDIR"
		$script2 50pourc $inputDIR $outputDIR

		echo -e "\t\t$script2 ALL $inputDIR $outputDIR"
		$script2 ALL $inputDIR $outputDIR

	done
done


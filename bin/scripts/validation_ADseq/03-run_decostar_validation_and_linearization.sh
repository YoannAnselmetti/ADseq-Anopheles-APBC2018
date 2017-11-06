#! /bin/bash

decostar="bin/software_libraries/DeCoSTAR/bin/DeCoSTAR"
paramDIR="data/validation_ADseq/DeCoSTAR/BLASTn"

READ="50pourc ALL"
ASSEMBLY="blastn"
SPECIES="Aalb-Anopheles_albimanus Aara-Anopheles_arabiensis Adir-Anopheles_dirus"

scriptDIR="bin/scripts/validation_ADseq/decostar/"

for tag in ASSEMBLY; do
	for spe in SPECIES; do
		speTAG=$(echo $spe|cut -d- -f1)
		speNAME=$(echo $spe|cut -d- -f2)
		for samp in READ; do
			$decostar parameter.file=$paramDIR"/"$speNAME"/"$samp"/INPUT_DeCoSTAR/decostar//validation_ADseq_"$speTAg"_"$samp"_Xtopo_spi20.param_"$tag".txt"
			$decostar parameter.file=$paramDIR"/"$speNAME"/"$samp"/INPUT_DeCoSTAR/decostar//validation_ADseq_"$speTAg"_"$samp"_Xtopo_ARt-DeCo_spi20.param_"$tag".txt"
		done
		$scriptDIR/run_postdecostar_validation_$speTAG.sh spi_20 $tag
	done
done
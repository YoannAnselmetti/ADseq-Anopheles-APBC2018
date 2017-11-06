#!/bin/bash
###
###   Goal:
###      Script to run all scripts for Validation of ADseq algo with BLASTn align
###
###
###   OUTPUT:
###      - compute stats to compare adjacencies prediction between Adseq, ARt-DeCo and BESST
###
###   Name: run_all_exp.sh                          Author: Yoann Anselmetti
###   Creation date: 2017/01/03                     Last modification: 2017/10/25
###

out="results/validation_ADseq/spi_20/stats/"
SCRIPTdir="bin/scripts/validation_ADseq/stats/01-compute_stats"

boolORI="with_ori without_ori"
assembly="blastn"
SPIlist="spi_20"

for ori in $boolORI; do
	for ass in $assembly; do
		for spi in $SPIlist; do
			outputDIR=$out/$ass/$ori
			mkdir -p $outputDIR
			echo $outputDIR

			echo "For species Aalb:"
			$SCRIPTdir/run_for_1exp.sh Anopheles_albimanus Aalb 50pourc $outputDIR data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/50pourc/Anopheles_albimanus/BESST_output/pass2/info-pass2.agp $spi $ass $ori
			$SCRIPTdir/run_for_1exp.sh Anopheles_albimanus Aalb ALL $outputDIR data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/ALL/Anopheles_albimanus/BESST_output/pass3/info-pass3.agp $spi $ass $ori
			echo "For species Aara:"
			$SCRIPTdir/run_for_1exp.sh Anopheles_arabiensis Aara 50pourc $outputDIR data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/50pourc/Anopheles_arabiensis/BESST_output/pass2/info-pass2.agp $spi $ass $ori
			$SCRIPTdir/run_for_1exp.sh Anopheles_arabiensis Aara ALL $outputDIR data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/ALL/Anopheles_arabiensis/BESST_output/pass2/info-pass2.agp $spi $ass $ori
			echo "For species Adir:"
			$SCRIPTdir/run_for_1exp.sh Anopheles_dirus Adir 50pourc $outputDIR data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/50pourc/Anopheles_dirus/BESST_output/pass5/info-pass5.agp $spi $ass $ori
			$SCRIPTdir/run_for_1exp.sh Anopheles_dirus Adir ALL $outputDIR data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/TRIMMOMATIC3/minia/ALL/Anopheles_dirus/BESST_output/pass5/info-pass5.agp $spi $ass $ori
		done
	done
done
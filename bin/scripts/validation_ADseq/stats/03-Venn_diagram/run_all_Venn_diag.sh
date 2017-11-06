#!/bin/bash
###
###   Goal:
###      Run all experiments for the Venn diagram
###   INPUT:

###   OUTPUT:
###      - Venn diagrams
###
###   Name: run_all_Venn_diag.sh     Author: Yoann Anselmetti
###   Creation date: 2017/02/19      Last modification: 2017/02/25
###




listassemb="blastn"
listspe="Aalb Aara Adir"
listseq="50pourc ALL"
listadj="TP CFP FP-CFP ALL FN"
listsupp="0.1 0.5 0.8"
scriptDIR="./bin/scripts/validation_ADseq/stats/03-Venn_diagram"

# listspe="Adir"
# listseq="ALL"
# listadj="CFP"
# listsupp="0.1"



for assemb in $listassemb; do
	vennWithoutORI="results/validation_ADseq/spi_20/stats/"$assemb"/without_ori/Venn_diagram"
	vennWithORI="results/validation_ADseq/spi_20/stats/"$assemb"/with_ori/Venn_diagram"
	outputWithORI="figures/Venn_diagram/"$assemb"/with_ori/"
	outputWithoutORI="figures/Venn_diagram/"$assemb"/without_ori/"
	for spe in $listspe; do
		for seq in $listseq; do
			for adj in $listadj; do
				for supp in $listsupp; do
					$scriptDIR/Venn_diagram_3comp.py $vennWithoutORI $spe $seq $adj $supp f $outputWithoutORI;
					$scriptDIR/Venn_diagram_3comp.py $vennWithORI $spe $seq $adj $supp t $outputWithORI;
				done
			done
		done
	done
done
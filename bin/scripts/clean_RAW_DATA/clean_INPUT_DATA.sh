#!/bin/bash
###
###   Goal:
###      Clean and standardize INPUT DATA available for the 18 Anopheles dataset
###
###   OUTPUT:
###      - Cleanned and standardized files for DeCoSTAR pipeline
###
###   Name: clean_IPUT_DATA.sh          Author: Yoann Anselmetti
###   Creation date: 2017/10/26         Last modification: 2017/10/26
###

scriptDIR="bin/scripts/clean_RAW_DATA/code"

nameID="data/INPUT_DATA/name_geneID_18Anopheles"
inputGT="data/INPUT_DATA/OG_CDS_newtrees"
intermGT="data/INPUT_DATA/unrooted_trees_nofilt.nwk"
outputGT="data/INPUT_DATA/unrooted_raw_trees.nwk"

FASTA="data/INPUT_DATA/FASTA/CDS"
inputGFF="data/INPUT_DATA/ORIGINAL_GFF"
outputGFF="data/INPUT_DATA/GFF"

if [ -d "$inputGT" ];
then
	echo ""$scriptDIR/clean_gene_trees.sh $inputGT $intermGT""
	$scriptDIR/clean_gene_trees.sh $inputGT $intermGT

	echo -e "\n"$scriptDIR/filter_species_in_gene_trees.py $intermGT $nameID $outputGT 1""
	$scriptDIR/filter_species_in_gene_trees.py $intermGT $nameID $outputGT 1

	echo -e "\n"$scriptDIR/transform_GFF_files.sh $FASTA $inputGFF $outputGFF""
	$scriptDIR/transform_GFF_files.sh $FASTA $inputGFF $outputGFF
else
	echo "Check if directory \"data/INPUT_DATA/OG_CDS_newtrees\" is uncompressed."
fi
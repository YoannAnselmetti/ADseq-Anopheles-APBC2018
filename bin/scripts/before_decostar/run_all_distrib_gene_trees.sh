#!/bin/bash

script="./bin/scripts/before_decostar/code/distrib_gene_trees.py"
output="figures/distrib_gene_trees/"

GTone="data/INPUT_DATA/unrooted_raw_trees.nwk"
GTtwo="data/GENE_TREES/unrooted_trees_filtered.nwk"
GTthree="data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk"

nameID="data/INPUT_DATA/name_geneID_18Anopheles"

Xmax=40



$script $GTone $nameID $output"01-no_filter" $Xmax

$script $GTtwo $nameID $output"02-filter" $Xmax

$script $GTthree $nameID $output"03-pNJ" $Xmax
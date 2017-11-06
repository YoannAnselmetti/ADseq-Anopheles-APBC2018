#!/bin/bash

script="bin/scripts/before_decostar/code/plot_hist_scaff_BESST.py"
output="figures/besst_score/"

SCAFFone="data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3"
SCAFFtwo="data/data_DeCoSTAR/scaff_BESST_DeCoSTAR"

TAGone="RAW"
TAGtwo="DeCoSTAR"

echo ""$script $SCAFFone $output$TAGone $TAGone""
$script $SCAFFone $output$TAGone $TAGone

echo ""$script $SCAFFtwo $output$TAGtwo $TAGtwo""
$script $SCAFFtwo $output$TAGtwo $TAGtwo

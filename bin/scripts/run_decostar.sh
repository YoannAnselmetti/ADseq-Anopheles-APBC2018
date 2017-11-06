#! /bin/bash

decostar="bin/software_libraries/DeCoSTAR/bin/DeCoSTAR"
paramDIR="data/data_DeCoSTAR/decostar"

$decostar parameter.file=$paramDIR"/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff.param.txt"
$decostar parameter.file=$paramDIR"/Xtopo_pNJ/DeCoSTAR_Anopheles_Xtopo-scaff.param.txt"
$decostar parameter.file=$paramDIR"/Xtopo_pNJ/DeCoSTAR_Anopheles_Xtopo+scaff.param.txt"
$decostar parameter.file=$paramDIR"/Xtopo_RAW/RAW_Anopheles.param.txt"
#!/bin/bash
#################################################################################
#################################################################################
###                                                                           ###
###   Goal:                                                                   ###
###      Build genome reference with gmap_build to align genes on PacBio CTG  ###
###                                                                           ###
###   INPUT:                                                                  ###
###      1- Gene ID                                                           ###
###         (AARA)                                                            ###
###      2- ORI file (annotFile or scaffFile)                                 ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/annotG_Aalb_ORI_filt) ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/scaff_Aalb_ORI_3) ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/CTG_Aalb_ORI) ###
###      3- MINIA file (annotFile or scaffFile)                               ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/annotG_Aalb_50pourc_filt) ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/scaff_Aalb_50pourc_3) ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/CTG_Aalb_50pourc_filt) ###
###      4- OUTPUT file                                                       ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/annotFile_Aalb_50pourc_allspecies) ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/scaff_Aalb_50pourc_3_allspecies) ###
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/CTG_Aalb_50pourc_allspecies) ###
###                                                                           ###
###   OUTPUT:	(RUN in ~5-10min => !!! Use 5Gb memory (during few sec) !!!)  ###
###      - Genome reference database for given species from PacBio scaffolds  ###
###                                                                           ###
###   Name: 07-fusion_MINIA_ORI_file_BLASTn.sh  Author: Yoann Anselmetti      ###
###   Creation date: 2016/04/19                 Last modification: 2016/11/17 ###
###                                                                           ###
#################################################################################
#################################################################################

buf=$4"_buf"
output=$(dirname $4)

mkdir -p $output

grep -v $1 $2 > $buf;
cat $3 $buf > $4
rm $buf

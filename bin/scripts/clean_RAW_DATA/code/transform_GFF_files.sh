#!/bin/bash
###
###   Goal:
###      Script to create association file between species name and ID
###      that will be use to analyze genome annotation and compute stats
###
###   INPUT:
###      1- Directory of FASTA files CDS (From WaterHouse website)
###         (data/INPUT_DATA/FASTA/CDS)
###      2- Directory containing INPUT GFF files
###         (data/INPUT_DATA/ORIGINAL_GFF)
###      3- Directory containing OUTPUT GFF files
###         (data/INPUT_DATA/GFF)
###
###   OUTPUT:
###      - Transform INPUT GFF files in the standard GFF3 format
###
###   Name: transform_GFF_files.sh        Author: Yoann Anselmetti
###   Creation date: 2015/10/28           Last modification: 2017/10/26
###                                                                          

mkdir -p $3

declare -A ID_NAME
for file in $(ls $1);
do
	ID=$(echo $file|cut -d. -f1)
	name=$(echo $file|cut -d. -f2)
	ID_NAME+=([$ID]=$name)
done 


# Course all GFF files unsorted in INPUT directory
for GFFfile in $(ls $2)
do
	ID=$(echo $GFFfile|cut -d. -f1)
	name=${ID_NAME["$ID"]}
	if [[ ! -z $name ]];
	then
		echo -e "\t=> "$name
		awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$11}' $2/$GFFfile > $3/$name".gff3"
	fi
done
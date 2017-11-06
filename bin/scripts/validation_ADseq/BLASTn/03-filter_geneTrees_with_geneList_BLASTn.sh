#!/bin/bash
###
###   Goal:
###      Filter gene trees containing gene present in file of genes to remove
###			=> 6 RUNS: 2: (50pourc & ALL)  x3: Aalb, Aara & Adir
###
###   INPUT:
###      1- File containing list of genes to remove
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/unmapG_Aalb_50pourc)
###      2- INPUT Gene trees file
###         (DECOSTAR/DATA_DeCoSTAR/GENE_TREES/ROOTED/rooted_trees_DeCoSTAR_Xtopo.nwk)
###      3- OUTPUT Gene tres file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/rooted_trees_DeCoSTAR_Xtopo_filt.nwk)
###
###   OUTPUT:
###      - Filtered gene trees file
###
###   Name: 03-filter_geneTrees_with_geneList_BLASTn.sh		Author: Yoann Anselmetti
###   Creation date: 2015/10/19								Last modification: 2016/11/17
###

geneNbTOT=$(wc -l $1)

mkdir -p $(dirname $3)

cp $2 $3

buf=$3"_buf"


GTnumber=$(wc -l $3 | cut -d" " -f1)
echo -e "\t=> "$GTnumber" gene trees before filtering"

i=1
while read line; do
	# echo "Processing gene "$i"/"$geneNbTOT;
	gene=$(echo $line | cut -d"|" -f1);
	grep -v $gene $3 > $buf;
	mv $buf $3;
	i=$(($i+1));
done < $1

GTnumber=$(wc -l $3 | cut -d" " -f1)
echo -e "\t=> "$GTnumber" gene trees after filtering\n"
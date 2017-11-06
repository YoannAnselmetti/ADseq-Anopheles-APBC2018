#!/bin/bash
###
###   Goal:
###      Concatenate raw gene trees files in 1 file, and clean it (remove empty gene trees and add "\n" where it is necessary).
###   INPUT:
###      1- Directory with unrooted gene trees
###         (data/INPUT_DATA/OG_CDS_newtrees)
###      2- OUTPUT gene trees file
###         (data/INPUT_DATA/unrooted_trees_nofilt.nwk)
###
###   OUTPUT:
###      - Gene trees file cleanned
###
###   Name: clean_gene_trees.sh                Author: Yoann Anselmetti
###   Creation date: 2015/09/22                Last modification: 2017/10/26
###

# After uncompressing "trees_ALL.tar.gz" obtained from WaterHouse website in DATA/DATA_WaterHouse/GENE_TREES/.
# We obtain the directory "OG_CDS_newtrees" with all unrooted gene trees obtained with RAxML in Neafsey et al., 2015, Science.

# CLEAN GENE TREES (Remove empty gene trees and add )
echo -n "Delete EMPTY gene trees and add \""
echo -n "\n\" between 2 gene trees if necessary (1 gene tree/line) ..."
# cat $1/* | sed 's/();//g' | sed 's/;(/;\n(/g' > $2
cat $1/* | sed 's/();//g' | sed 's/;(/;\n(/g' > $2
echo "DONE"
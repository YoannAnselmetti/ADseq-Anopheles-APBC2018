# Cedric Chauve, December 20, 2016
# Linearizing and analyzing evolution for the genomes of the decostar datasets

#! /bin/bash

scriptLINEAR=bin/scripts/post_decostar/code/linearize_assign_scj.sh
resultsDIR=results/decostar

### For DeClone support threshold to 0.1 for linearization
seuil=0.1

GENEfile=data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF

### WGtopo+scaff

DIR=$resultsDIR/WGtopo+scaff
longPREF=$DIR/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1
PREF=DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_$seuil

$scriptLINEAR ${longPREF}.reconciliations.newick \
            ${longPREF}.speciesTree.newick \
            ${longPREF}.adjacencies.txt \
            $GENEfile \
            $seuil \
            $DIR/ \
            ${PREF} \
            M1 \
            python



### Xtopo+scaff
DIR=$resultsDIR/Xtopo+scaff
longPREF=$DIR/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1
PREF=DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_$seuil

$scriptLINEAR ${longPREF}.reconciliations.newick \
            ${longPREF}.speciesTree.newick \
            ${longPREF}.adjacencies.txt \
            $GENEfile \
            $seuil \
            $DIR \
            ${PREF} \
            M1 \
            python



### Xtopo-scaff
DIR=$resultsDIR/Xtopo-scaff
longPREF=$DIR/DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1
PREF=DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1_$seuil

$scriptLINEAR ${longPREF}.reconciliations.newick \
            ${longPREF}.speciesTree.newick \
            ${longPREF}.adjacencies.txt \
            $GENEfile \
            $seuil \
            $DIR \
            ${PREF} \
            M1 \
            python



### Xtopo_RAW
DIR=$resultsDIR/Xtopo_RAW
longPREF=$DIR/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1
PREF=DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_$seuil

$scriptLINEAR ${longPREF}.reconciliations.newick \
            ${longPREF}.speciesTree.newick \
            ${longPREF}.adjacencies.txt \
            $GENEfile \
            $seuil \
            $DIR \
            ${PREF} \
            M1 \
            python
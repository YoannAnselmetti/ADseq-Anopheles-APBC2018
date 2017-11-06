#!/bin/bash
###
###   Goal:
###      Script to run all scripts for Validation of ADseq algo with BLASTn align
###
###   INPUT:
###      1- Species name
###         (Anopheles_albimanus / Anopheles_arabiensis / Anopheles_dirus)
###      2- Species TAG
###         (Aalb / Aara / Adir)
###      3- Reads sampling TAG_ID
###         (50pourc / ALL)
###      4- ORIGINAL/REFERENCE genome assembly
###         (data/INPUT_DATA/FASTA/SCAFF/Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa)
###         (data/INPUT_DATA/FASTA/SCAFF/Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa)
###         (data/INPUT_DATA/FASTA/SCAFF/Anopheles-dirus-WRAIR2_SCAFFOLDS_AdirW1.fa)
###      5- MINIA genome assembly
###        (data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_albimanus/50pourc/minia_k75_m3_Aalb_50pourc.contigs.fa)
###		   (data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_albimanus/ALL/minia_k83_m3_Aalb_ALL.contigs.fa)
###        (minia_k59_m3_Aara_50pourc.contigs.fa / minia_k72_m3_Aara_ALL.contigs.fa)
###        (minia_k63_m3_Adir_50pourc.contigs.fa / minia_k75_m3_Adir_ALL.contigs.fa)
###      6- Number of the script where pipeline have to begin(Avoid to run all scripts if the X first are done)
###        (3 => Begin to filter of gene trees with BLASTn results)
###      7- GENE species TAG ID
###         (AALB / AARA / ADIR)
###      8- assembly TAG
###         (minia / blastn)
###
###   OUTPUT:
###      - BLASTn results of CTG from MINIA assembly mapped on CTG from REF assembly to the gene ID
###
###   Name: 00-run_all_scripts_BLASTn.sh             Author: Yoann Anselmetti
###   Creation date: 2016/11/17                      Last modification: 2017/10/26
###

##################
### INPUT DATA ###
##################

# Reference genome assemblies
FASTAdir="data/INPUT_DATA/FASTA/SCAFF"
# Directories containing sequencing data (scaffolding adjacencies) obtained with Bowtie2+BESST
SCAFFdir=""
if [ $8 = "minia" ]; then
	SCAFFdir="data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_k50/TRIMMOMATIC3/minia"
elif [ $8 = "blastn" ]; then
	SCAFFdir="data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_k50/TRIMMOMATIC3/blastn"
fi
SCAFFdirORI="data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_ALL/TRIMMOMATIC3/ALL"
# Annotation gene file obtained from reference genome assemblies data
GENEfile="data/data_DeCoSTAR/GENE_file"
# Gene trees obtained from INPUT data
INPUTtrees="data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk"
RAWtrees="data/INPUT_DATA/OG_CDS_newtrees"
# Path to scripts repertory
SCRIPTdir="bin/scripts/validation_ADseq/BLASTn"


###################
### OUTPUT DATA ###
###################

# Path for OUTPUT directories
ADSEQdir="data/validation_ADseq/DeCoSTAR/BLASTn/"$1"/"$3
BLASTNdir="data/validation_ADseq/BLASTn"
BLASTNalign=$BLASTNdir"/"$1"/"$3"/align_"$2"_ORI_"$3"_allparams.tab"


# Files of species used for simulation of genome fragmentation for validation of ADseq
annotGfile=$ADSEQdir"/annotG_"$2"_"$3
annotGfileFILT=$ADSEQdir"/annotG_"$2"_"$3"_filt"
CTGfile=$ADSEQdir"/CTG_"$2"_"$3
CTGfileFILT=$ADSEQdir"/CTG_"$2"_"$3"_filt"
assoCTGfile=$ADSEQdir"/assoCTG_"$2"_"$3
unmapGfile=$ADSEQdir"/unmapG_"$2"_"$3
preSCAFFfile=$ADSEQdir"/preScaff_"$2"_"$3"_3_"$8
preSCAFFfileMODIF=$ADSEQdir"/preScaff_"$2"_"$3"_3_"$8"_modif"
SCAFFfile=$ADSEQdir"/scaff_"$2"_"$3"_3_"$8

# Files for genome assemblies for others species and comparison with genome assembly of species with genome fragmentation simulation
annotGfileORI=$ADSEQdir"/annotG_"$2"_ORI_filt"
CTGfileORI=$ADSEQdir"/CTG_"$2"_ORI"
preSCAFFfileORI=$ADSEQdir"/preScaff_ORI_3"
SCAFFfileORI=$ADSEQdir"/scaff_"$2"_ORI_3"
spe=$(echo $1 |sed 's/_/-/g')
newSCAFFfile="data/validation_ADseq/FASTA/SCAFF/blastn/"$3"/"$1"/"$spe"-blastn_"$2"_"$3".fa"

# Files to execute DeCo* with ADseq algo
INPUTdecostar=$ADSEQdir"/INPUT_DeCoSTAR"
mkdir -p $INPUTdecostar
OUTPUTtrees=$INPUTdecostar"/rooted_trees_DeCoSTAR_filt.nwk"
annotGfileFINAL=$INPUTdecostar"/annotG_"$2"_"$3"_allspecies"
CTGfileFINAL=$INPUTdecostar"/CTG_"$2"_"$3"_allspecies"
SCAFFfileFINAL=$INPUTdecostar"/scaff_"$2"_"$3"_3_allspecies_"$8

# DeCoSTAR INPUT files
DECOSTARdir=$INPUTdecostar"/decostar"
ADJfileADseq=$DECOSTARdir"/adjacencies_"$2"_"$3"_"$8
ADJfileAD=$DECOSTARdir"/adjacencies_"$2"_"$3"_ARt-DeCo"
TREESdir=$DECOSTARdir"/DeCoSTAR_Anopheles_Xtopo_gene_trees"
PATHtree=$DECOSTARdir"/DeCoSTAR_Anopheles_Xtopo_gene_trees"
TREESdistrib=$DECOSTARdir"/distrib_DeCoSTAR_Anopheles_Xtopo_gene_trees.txt"


###############################################################################################################################
### SCRIPTS for ADseq algo validation by simulation of genome fragmentation with genome annotation made by BLASTn alignment ###
###############################################################################################################################

# Script 1: Align MINIA contigs on reference genome assembly to annotate MINIA contigs 
if [ 1 -ge $6 ]; then
	echo "$SCRIPTdir/01-run_blastn.sh $4 $5 $BLASTNalign"
	$SCRIPTdir/01-run_blastn.sh $4 $5 $BLASTNalign 
fi



# Script 2: filter alignements BLASTn results to create genome annotation
if [ 2 -ge $6 ]; then
	echo "$SCRIPTdir/02-filter_BLASTn_results.py $GENEfile $BLASTNalign $CTGfile $annotGfile $assoCTGfile $unmapGfile $1 90 90 f $5 $newSCAFFfile $4"
	$SCRIPTdir/02-filter_BLASTn_results.py $GENEfile $BLASTNalign $CTGfile $annotGfile $assoCTGfile $unmapGfile $1 90 90 f $5 $newSCAFFfile $4
fi



if [ 3 -ge $6 ]; then
	echo "$SCRIPTdir/03-filter_geneTrees_with_geneList_BLASTn.sh $unmapGfile $INPUTtrees $OUTPUTtrees"
	$SCRIPTdir/03-filter_geneTrees_with_geneList_BLASTn.sh $unmapGfile $INPUTtrees $OUTPUTtrees
fi



if [ 4 -ge $6 ]; then
	echo "$SCRIPTdir/04-filter_Gene_with_geneTrees_BLASTn.py $annotGfile $OUTPUTtrees $annotGfileFILT @ prefix $1"
	$SCRIPTdir/04-filter_Gene_with_geneTrees_BLASTn.py $annotGfile $OUTPUTtrees $annotGfileFILT @ prefix $1
	echo "$SCRIPTdir/04-filter_Gene_with_geneTrees_BLASTn.py $GENEfile $OUTPUTtrees $annotGfileORI @ prefix $1"
	$SCRIPTdir/04-filter_Gene_with_geneTrees_BLASTn.py $GENEfile $OUTPUTtrees $annotGfileORI @ prefix $1
fi



if [ 5 -ge $6 ]; then
	echo "$SCRIPTdir/05a-create_CTG_file_BLASTn.py $annotGfileFILT $CTGfileFILT $1"
	$SCRIPTdir/05a-create_CTG_file_BLASTn.py $annotGfileFILT $CTGfileFILT $1
	echo "$SCRIPTdir/05b-create_CTG_file_ORI_BLASTn.py $annotGfileORI $FASTAdir $CTGfileORI"
	$SCRIPTdir/05b-create_CTG_file_ORI_BLASTn.py $annotGfileORI $FASTAdir $CTGfileORI
fi



if [ 6 -ge $6 ]; then

	# For fragmentated genome
	echo "$SCRIPTdir/06a-create_scaff_adj_prefile_BLASTn.py $SCAFFdir"/"$3 1000000000 3 $preSCAFFfile $1"
	$SCRIPTdir/06a-create_scaff_adj_prefile_BLASTn.py $SCAFFdir"/"$3 1000000000 3 $preSCAFFfile $1
	if [ $8 = "minia" ]; then
		# Devrait disparaitre une fois que le pipeline avec fusion complète des contigs sera fait.
		echo "$SCRIPTdir/06b-modif_scaff_adj_prefile_BLASTn.py $preSCAFFfile $assoCTGfile $preSCAFFfileMODIF"
		$SCRIPTdir/06b-modif_scaff_adj_prefile_BLASTn.py $preSCAFFfile $assoCTGfile $preSCAFFfileMODIF
		echo "$SCRIPTdir/06b-create_scaff_adj_file_final_BLASTn.py $preSCAFFfileMODIF $CTGfileFILT $SCAFFfile"
		$SCRIPTdir/06b-create_scaff_adj_file_final_BLASTn.py $preSCAFFfileMODIF $CTGfileFILT $SCAFFfile
	elif [ $8 = "blastn" ]; then
		echo "$SCRIPTdir/06b-create_scaff_adj_file_final_BLASTn.py $preSCAFFfile $CTGfileFILT $SCAFFfile"
		$SCRIPTdir/06b-create_scaff_adj_file_final_BLASTn.py $preSCAFFfile $CTGfileFILT $SCAFFfile
	fi

	# For genomes not fragmentated (ORI)
	echo "$SCRIPTdir/06a-create_scaff_adj_prefile_BLASTn.py $SCAFFdirORI 1000000000 3 $preSCAFFfileORI ALL"
	$SCRIPTdir/06a-create_scaff_adj_prefile_BLASTn.py $SCAFFdirORI 1000000000 3 $preSCAFFfileORI ALL
	echo "$SCRIPTdir/06b-create_scaff_adj_file_final_BLASTn.py $preSCAFFfileORI $CTGfileORI $SCAFFfileORI"
	$SCRIPTdir/06b-create_scaff_adj_file_final_BLASTn.py $preSCAFFfileORI $CTGfileORI $SCAFFfileORI
fi


# Produce data files necessary for DeCoSTAR (merge data of species with fragmented genome simulation and the remaining species of the dataset)
if [ 7 -ge $6 ]; then
	echo "$SCRIPTdir/07-fusion_MINIA_ORI_file_BLASTn.sh $7 $annotGfileORI $annotGfileFILT $annotGfileFINAL"
	$SCRIPTdir/07-fusion_MINIA_ORI_file_BLASTn.sh $7 $annotGfileORI $annotGfileFILT $annotGfileFINAL
	echo "$SCRIPTdir/07-fusion_MINIA_ORI_file_BLASTn.sh $7 $CTGfileORI $CTGfileFILT $CTGfileFINAL"
	$SCRIPTdir/07-fusion_MINIA_ORI_file_BLASTn.sh $7 $CTGfileORI $CTGfileFILT $CTGfileFINAL
	echo "$SCRIPTdir/07-fusion_MINIA_ORI_file_BLASTn.sh $7 $SCAFFfileORI $SCAFFfile $SCAFFfileFINAL"
	$SCRIPTdir/07-fusion_MINIA_ORI_file_BLASTn.sh $7 $SCAFFfileORI $SCAFFfile $SCAFFfileFINAL
fi


# Create adjacencies file used as INPUT of DeCoSTAR
if [ 8 -ge $6 ]; then
	echo "$SCRIPTdir/08-create_ADJfile_for_DeCoSTAR_BLASTn.py $annotGfileFINAL $SCAFFfileFINAL $ADJfileADseq $ADJfileAD @ prefix"
	$SCRIPTdir/08-create_ADJfile_for_DeCoSTAR_BLASTn.py $annotGfileFINAL $SCAFFfileFINAL $ADJfileADseq $ADJfileAD @ prefix
fi


# Create gene trees directory for DeCoSTAR
if [ 9 -ge $6 ]; then
	echo "$SCRIPTdir/09-write_1tree_per_file_for_DeCoSTAR_BLASTn.py $OUTPUTtrees $RAWtrees $TREESdir $PATHtree $TREESdistrib @ prefix"
	$SCRIPTdir/09-write_1tree_per_file_for_DeCoSTAR_BLASTn.py $OUTPUTtrees $RAWtrees $TREESdir $PATHtree $TREESdistrib @ prefix
fi
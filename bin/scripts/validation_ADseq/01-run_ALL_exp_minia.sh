#!/bin/bash
###
###   Goal:
###      Script to run all scripts to produce a new genome assembly with minia
###      for the 3 species: An. albimanus, An. arabiensis and An. dirus
###      and the two reads sampling: 100% (ALL) and 50% (50pourc)
###
###
###   Name: 01-run_ALL_exp_minia.sh             Author: Yoann Anselmetti
###   Creation date: 2016/11/17                 Last modification: 2017/10/27
###


redFASTQ="bin/scripts/validation_ADseq/minia/reduce_pairedFASTQ_files.py"
calcN50="bin/scripts/validation_ADseq/minia/calcN50_fasta.pl"
sampling=0.5
SPECIES="Aalb-Anopheles_albimanus Aara-Anopheles_arabiensis Adir-Anopheles_dirus"

FASTQdir="data/DATA_SEQ/FASTQ/"

kmerG="bin/software_libraries/kmergenie-1.7016/kmergenie"
kmerDIR="data/validation_ADseq/kmergenie/"

minia="bin/software_libraries/minia-2.0.3-Source/build/minia"
miniaDIR="data/validation_ADseq/FASTA/SCAFF/minia/"

logALL=$kmerDIR/$speNAME/ALL/log_file.txt
logHALF=$kmerDIR/$speNAME/50pourc/log_file.txt

for spe in SPECIES;
do
	speTAG=$(echo $spe|cut -d- -f1)
	speNAME=$(echo $spe|cut -d- -f2)

	echo "STEP1: Sampling of FASTQ files of species "$speNAME":"
	$redFASTQ $sampling data/DATA_SEQ/FASTQ/ALL/$speNAME

	echo "STEP2: Get list of FASTQ files of species "$speNAME":"
	ls -1 $FASTQdir/ALL/$speNAME/SRX*/*/* > $FASTQdir/list_FASTQ/list_ALL_FASTQ_$speTAG".txt"
	ls -1 $FASTQdir/50pourc/$speNAME/SRX*/*/* > $FASTQdir/list_FASTQ/list_50pourc_FASTQ_$speTAG".txt"

	echo "STEP3: Compute value of the kmer size with KmerGenie for the assembly of "$speNAME" with minia assembly tool:"
	$kmerG $FASTQdir/list_FASTQ/list_ALL_FASTQ_$speTAG".txt" -o $kmerDIR/$speNAME/ALL/kmergenie_$speTAG_ALL_histograms > $logALL
	$kmerG $FASTQdir/list_FASTQ/list_50pourc_FASTQ_$speTAG".txt" -o $kmerDIR/$speNAME/50pourc/kmergenie_$speTAG_50pourc_histograms > $logHALF

	kmerALL=$(tail -1 $logALL |awk '{print $NF}')
	kmerHALF=$(tail -1 $logHALF |awk '{print $NF}')

	echo "STEP4: Genome assembly of "$speNAME" with minia tool:"
	mkdir -p data/DATA_SEQ/minia/Anopheles_albimanus/ALL;
	$minia -in $FASTQdir/list_FASTQ/list_ALL_FASTQ_$speTAG".txt" -kmer-size $kmerALL -abundance-min 3 -out $miniaDIR/$speNAME/ALL/minia_k$kmerALL_m3_$speTAG_ALL
	$minia -in $FASTQdir/list_FASTQ/list_50pourc_FASTQ_$speTAG".txt" -kmer-size $kmerHALF -abundance-min 3 -out $miniaDIR/$speNAME/ALL/minia_k$kmerHALF_m3_$speTAG_50pourc

	echo -e "\n\nAssembly statistics for species "$speNAME":\n"
	# ASSEMBLY STATISTICS
	
	# N50(bp) + (#CTG>=N50):
	echo -e "\nN50 statistics:\n"
	echo -e "\tALL:\n"
	$calcN50 -f $miniaDIR/$speNAME/ALL/minia_k$kmerALL_m3_$speTAG_ALL.contigs.fa -l 0
	echo -e "\t50pourc:\n"
	$calcN50 -f $miniaDIR/$speNAME/50pourc/minia_k$kmerHALF_m3_$speTAG_50pourc.contigs.fa -l 0
	
	# #CTG:
	# 	=> See in output file ".o" line with TAG "nb_contigs".
	echo -e "\nNumber of CTG:\n"
	echo -e "\tALL:\n"
	grep -c \> $miniaDIR/$speNAME/ALL/minia_k$kmerALL_m3_$speTAG_ALL.contigs.fa
	echo -e "\t50pourc:\n"
	grep -c \> $miniaDIR/$speNAME/50pourc/minia_k$kmerHALF_m3_$speTAG_50pourc.contigs.fa
	
	# #CTG>=1000bp:
	echo -e "\nNumber of CTG with size>1000bp:\n"
	echo -e "\tALL:\n"
	grep \> $miniaDIR/$speNAME/ALL/minia_k$kmerALL_m3_$speTAG_ALL.contigs.fa | awk -F"__" '{if ($3>=1000) sum+=1} END {print sum}'
	echo -e "\t50pourc:\n"
	grep \> $miniaDIR/$speNAME/50pourc/minia_k$kmerHALF_m3_$speTAG_50pourc.contigs.fa | awk -F"__" '{if ($3>=1000) sum+=1} END {print sum}'

	# Assembly size:
	echo -e "\nAssembly size:\n"
	echo -e "\tALL:\n"
	# 	=> See in output file ".o" line with TAG "nt_assembled"
	grep \> $miniaDIR/$speNAME/ALL/minia_k$kmerALL_m3_$speTAG_ALL.contigs.fa | awk -F"__" '{sum+=$3} END {print sum}'
	echo -e "\t50pourc:\n"
	grep \> $miniaDIR/$speNAME/50pourc/minia_k$kmerHALF_m3_$speTAG_50pourc.contigs.fa | awk -F"__" '{sum+=$3} END {print sum}'

done
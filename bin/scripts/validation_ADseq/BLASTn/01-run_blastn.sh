#!/bin/bash
###
###   Goal:
###      Apply BLASTn as query CTG of MINIA genomeassembly and Suject SCAFF of ORIGINAL/REFERENCE genome assembly
###      	=> RUN 6 experiments: 2 reads sampling (50% & ALL) x 3 species (Aalb, Aara & Adir)
###
###   INPUT:
###      1- ORIGINAL/REFERENCE genome assembly
###         (DECOSTAR/DATA_DeCoSTAR/INPUT_DATA/FASTA/SCAFF/Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa)
###      2- MINIA genome assembly
###         (Validation_ADseq/coverage/minia/Anopheles_arabiensis/50pourc/minia_k59_m3_Aara_50pourc.contigs.fa)
###      3- OUTPUT file results of BLASTn (megablast)
###         (Validation_ADseq/coverage/BLASTn/Anopheles_arabiensis/50pourc/align_Aara_ORI_50pourc_allparams.tab)
###
###   OUTPUT:
###      - BLASTn results of CTG from MINIA assembly mapped on CTG from REF assembly to the gene ID
###
###   Name: 01-run_blastn.sh                         Author: Yoann Anselmetti
###   Creation date: 2016/07/25                      Last modification: 2016/11/17
###

mkdir -p $(dirname $3)
BLASTn=bin/software_libraries/ncbi-blast-2.4.0+/bin/blastn
# TODO => Allow several format for BLASTn (BLAST.XML / BLAST default format / ...) 1E-10
$BLASTn -task megablast -subject $1 -query $2 -evalue 1E-10 -outfmt 6" qseqid qlen sseqid slen pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore" -out $3
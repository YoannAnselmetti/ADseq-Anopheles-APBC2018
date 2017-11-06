#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create contigs file with gene located to the extremities of contigs
###
###   INPUT:
###      1- Annotation gene file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/annotG_Aara_ORI_filt)
###      2- Genome Assembly FASTA files directory
###         (DECOSTAR/DATA_DeCoSTAR/INPUT_DATA/FASTA/SCAFF/)
###      3- Contigs OUTPUT file giving info on genes to extremities
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/CTG_Aara_ORI)
###
###   OUTPUT:
###      - Contig file with genes located to extremities
###
###   Name: 05b-create_CTG_file_ORI_BLASTn.py  Author: Yoann Anselmetti
###   Creation date: 2016/02/15                Last modification: 2017/03/15
###

from sys import argv, stdout
from re import search
from os import close, path, mkdir, listdir
from datetime import datetime
from collections import namedtuple   #New in version 2.6


def mkdir_p(dir_path):
	try:
		makedirs(dir_path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and path.isdir(dir_path):
			pass
		else:
			raise

################
###   MAIN   ###
################
if __name__ == '__main__':

	start_time = datetime.now()

	# Recovery of input parameters
	annotG_file=argv[1]
	SCAFF_dir=argv[2]
	CTG_file=argv[3]

	OUTPUT_DIR=path.dirname(path.realpath(CTG_file))
	# Create OUTPUT_DIR if not existing
	if not path.exists(OUTPUT_DIR):
		mkdir_p(OUTPUT_DIR)

	# To be sure than directory have no "/" to the end of the path
	SCAFF_dir=path.normpath(SCAFF_dir)

	CTG=namedtuple("CTG",["spe","size"])

#################################################################
### Store list of contig size/species in dict_spe_listCTGsize ###
#################################################################
	dict_CTG={}
	FASTA_list=listdir(SCAFF_dir)
	print "Storing CTG size of species:"
	# Browse list of Genome assemblies to FASTA file format
	for FASTA in sorted(FASTA_list):
		species=""
		r=search("^([^-]*)-([^-]*)-.*\.fa$",FASTA)
		if r:
			genus=r.group(1)
			spe=r.group(2)
			species=genus+"_"+spe
		else:
			exit("!!! ERROR, FASTA file name: "+FASTA+" is incorrectly written !!!")

		print "\t"+species
		# Browse current FASTA file to get list of contigs with its size 
		fasta_file=open(SCAFF_dir+"/"+FASTA)
		for line in fasta_file:
			r_line=search("^>([^ ]*) .*:([0-9]*):[^:\n]*\n$",line)
			if r_line:
				contig=r_line.group(1)
				size=r_line.group(2)
	            # print "\t\t"+contig+" - "+size+" bp"
				if contig in dict_CTG:
					exit("!!! ERROR, there are two contigs that have the same ID ("+contig+") !!!")

				ctg=CTG(species,size)
				dict_CTG[contig]=ctg

	fasta_file.close()

####################################################################
### BROWSE ANNOTATION GENE FILE TO WRITE GENE TO CTG EXTREMITIES ###
####################################################################
	print "\nBrowse gene annotation file to write CTG file...",
	dict_spe_newAdj={}
	input_file=open(annotG_file,'r')
	output_file=open(CTG_file,'w')
	output_file.write("#species\tctg\tctg_size\tctg_gene_nb\t5'_gene_family\t5'_gene\torientation_5'_gene\tstart_5'_gene\t3'_gene_family\t3'_gene\torientation_3'_gene\tend_3'_gene\n")
	gf1=""
	gene1=""
	ori1=""
	start1=""
	str_spe=""
	str_ctg=""
	str_gf=""
	str_gene=""
	str_ori="?"
	str_stop="?"
	size=""
	gene_nb=0
	for line in input_file:
		r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
		if r:
			cur_spe=r.group(1)
			cur_ctg=r.group(2)
			cur_gf=r.group(3)
			cur_gene=r.group(4)
			cur_ori=r.group(5)
			cur_start=r.group(6)
			cur_end=r.group(7)

			if cur_spe!="species":
				if str_spe=="":
					gf1=cur_gf
					gene1=cur_gene
					ori1=cur_ori
					start1=cur_start
				else:
					if (str_spe!=cur_spe) or (str_ctg!=cur_ctg):
						# print "str_ctg: "+str_ctg
						if str_ctg in dict_CTG:
						  output_file.write(str_spe+"\t"+str_ctg+"\t"+dict_CTG[str_ctg].size+"\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
						# For contig that are not present in SCAFF directory (Genome assemblies)
						else:
						  output_file.write(str_spe+"\t"+str_ctg+"\t?\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
						gene_nb=0
						gf1=cur_gf
						gene1=cur_gene
						ori1=cur_ori
						start1=cur_start
				str_spe=cur_spe
				str_ctg=cur_ctg
				str_gf=cur_gf
				str_gene=cur_gene
				str_ori=cur_ori
				str_stop=cur_end
				gene_nb+=1
	input_file.close()
	print "DONE"

	# Write last contig in CTG_file
	if str_ctg in dict_CTG:
	  output_file.write(str_spe+"\t"+str_ctg+"\t"+dict_CTG[str_ctg].size+"\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
	# For contig that are not present in SCAFF directory (Genome assemblies)
	else:
	  output_file.write(str_spe+"\t"+str_ctg+"\t?\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
	output_file.close()

	end_time = datetime.now()
	print('\nDuration: {}'.format(end_time - start_time))

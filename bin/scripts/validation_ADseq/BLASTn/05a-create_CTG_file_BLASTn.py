#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create contigs file with gene located to the extremities of contigs
###
###   INPUT:
###      1- Annotation gene file
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_arabiensis/50pourc/annotG_Aara_50pourc_filt)
###      2- Contigs OUTPUT file giving info on genes to extremities
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_arabiensis/50pourc/CTG_Aara_50pourc_filt)
###      3- Species name
###         (Anopheles_arabiensis)
###
###   OUTPUT:
###      - Contig file with genes located to extremities
###
###   Name: 05a-create_CTG_file_BLASTn.py      Author: Yoann Anselmetti
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
	CTG_file=argv[2]
	species=argv[3]

	OUTPUT_DIR=path.dirname(path.realpath(CTG_file))
	# Create OUTPUT_DIR if not existing
	if not path.exists(OUTPUT_DIR):
		mkdir_p(OUTPUT_DIR)

####################################################################
### BROWSE ANNOTATION GENE FILE TO WRITE GENE TO CTG EXTREMITIES ###
####################################################################
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
					print cur_spe
					gf1=cur_gf
					gene1=cur_gene
					ori1=cur_ori
					start1=cur_start
				else:
					if (str_spe!=cur_spe) or (str_ctg!=cur_ctg):
						size=str_ctg.split('__')[2]
						# print "str_ctg: "+str_ctg
						output_file.write(str_spe+"\t"+str_ctg+"\t"+size+"\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
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

	# Write last contig in CTG_file
	size=str_ctg.split('__')[2]
	output_file.write(str_spe+"\t"+str_ctg+"\t"+size+"\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
	output_file.close()

	end_time = datetime.now()
	print('\nDuration: {}'.format(end_time - start_time))

#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      parse and filter GENE file with gene trees file (and rewrite GENE family ID) to obtain a file with gene info for ARt-DeCo_seq instance creation
###         => 12 RUNS: 2 (ORI and spe)  x2: (50pourc & ALL)  x3: Aalb, Aara & Adir
###
###   INPUT:
###      1- INPUT annotation GENE file
###         (DECOSTAR/data_files_DeCoSTAR/GENE_file)
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/annotG_Aara_50pourc)
###      2- Gene trees file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/rooted_trees_DeCoSTAR_Xtopo_filt.nwk)
###      3- OUTPUT annotation gene file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/annotG_Aara_ORI_filt)
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/annotG_Aara_50pourc_filt)
###      4- Character separator between species name and gene ID
###         (@)
###      5- prefix/postfix boolean
###         (prefix or postfix)
###      6- Species processed
###         (Anôpheles_arabiensis)
###
###   OUTPUT:
###      - Create annotation gene file for ARt-DeCo_seq instance creation
###
###   Name: 04-filter_Gene_with_geneTrees_BLASTn.py  Author: Yoann Anselmetti
###   Creation date: 2015/11/11                      Last modification: 2016/11/17
###

from sys import argv
from re import search, match
from os import close, path, makedirs, listdir
from collections import namedtuple   #New in version 2.6
import subprocess
from datetime import datetime
from ete3 import Tree


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   GENE_file=argv[1]
   GT_file=argv[2]
   OUTPUT_annot_file=argv[3]
   separator=argv[4]
   order_bool=argv[5]
   select_spe=argv[6]

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_annot_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   GENE=namedtuple("GENE",["spe","ctg","gf","id","ori","start","stop"])


############################################
### INDEXATION OF GENES INFOS BY GENE ID ###
############################################
   print "INDEX of gene info by gene ID...",
   dict_ID_gene={}
   with open(GENE_file,'r') as ortho_file:
      for line in ortho_file:
         # If INPUT file: DECOSTAR/data_files_DeCoSTAR/GENE_file
         r_decostar=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         # Format of annotation gene file obtained with BLASTn
         r_blast=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         species=""
         contig=""
         GF_ID=""
         gene_ID=""
         gene_ori=""
         start_pos=""
         end_pos=""
         if r_decostar:
            species=r_decostar.group(1)
            contig=r_decostar.group(2)
            GF_ID=r_decostar.group(3)
            gene_ID=r_decostar.group(4)
            gene_ori=r_decostar.group(5)
            start_pos=r_decostar.group(6)
            end_pos=r_decostar.group(7)
            exon_nb=r_decostar.group(8)
            exon_list=r_decostar.group(9)
         elif r_blast:
            species=r_blast.group(1)
            contig=r_blast.group(2)
            GF_ID=r_blast.group(3)
            gene_ID=r_blast.group(4)
            gene_ori=r_blast.group(5)
            start_pos=r_blast.group(6)
            end_pos=r_blast.group(7)
         else:
            exit("\n!!! ERROR, in line:\n\t"+line+" of file "+GENE_file+" is not to the expected format\n!!!")

         if species!="species":
            # print GENE(species,contig,gene_ID,start_pos)
            gene_info=GENE(species,contig,GF_ID,gene_ID,gene_ori,start_pos,end_pos)
            dict_ID_gene[gene_ID]=gene_info
            # print gene_ID
   ortho_file.close()
   print "DONE"



#############################################################################################################
### BROWSE GENE TREES FILE AND GET GENE TO FILTER GENE IN ANNOTATION GENE FILE AND MODIFIY GENE FAMILY ID ###
#############################################################################################################
   print "Browse gene trees file "+GT_file+" to filter gene annotation file and add Gene Family ID...",
   list_genes=list()
   gene_remaining_all=0
   gene_remaining_spe=0
   # Read gene trees file to filter genes that are not present in annotation gene file
   input_trees=open(GT_file,"r")
   output_annot=open(OUTPUT_annot_file,"w")
   for tree_line in input_trees:
      tree=Tree(tree_line)
      # Get list of extant genes in current gene tree
      for leaf in tree.get_leaf_names():
         # print leaf
         gene=""
         spe=""
         if separator in leaf:
            if order_bool=="prefix":
               gene=leaf.split(separator)[1]
               spe=leaf.split(separator)[0]
            elif order_bool=="postfix":
               gene=leaf.split(separator)[0]
               spe=leaf.split(separator)[1]
            else:
               exit("ERROR, parameter 7 should be equal to \"postfix\" or \"prefix\" !!!")
         else:
            gene=leaf
         # print gene

         if gene in dict_ID_gene:
            if select_spe==spe:
               gene_remaining_spe+=1
            gene_remaining_all+=1
            list_genes.append(gene)
            output_annot.write(dict_ID_gene[gene].spe+"\t"+dict_ID_gene[gene].ctg+"\t"+dict_ID_gene[gene].gf+"\t"+dict_ID_gene[gene].id+"\t"+dict_ID_gene[gene].ori+"\t"+dict_ID_gene[gene].start+"\t"+dict_ID_gene[gene].stop+"\n")
            dict_ID_gene.pop(gene,None)

   print "DONE"
   print "\t=> "+str(len(dict_ID_gene))+" genes have been removed from GENE file "+GENE_file+" cause not present in gene trees file "+GT_file+" !!!"
   print "\t=> "+str(gene_remaining_spe)+" gene number remaining in species "+select_spe
   print "\t=> "+str(gene_remaining_all)+" gene number remaining in the dataset"
   output_annot.close()
   input_trees.close()

   # Write 1st line legend of annotation gene file
   output_buffer=open("buffer_file","w")
   output_buffer.write("species\tctg\tgene_family\tgene\tori_gene\tstart_gene\tend_gene\n")
   output_buffer.close()

   # Sort OUTPUT annotation gene file by species, then contigs and genes.
   command_line="sort -k1d,1d -k2d,2d -k6n,6n "+OUTPUT_annot_file+" >> buffer_file; mv buffer_file "+OUTPUT_annot_file
   subprocess.call(command_line,shell=True)

   # Get and print execution time of the script
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))

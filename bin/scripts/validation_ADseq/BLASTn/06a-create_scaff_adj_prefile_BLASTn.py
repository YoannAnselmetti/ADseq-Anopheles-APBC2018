#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###       Create scaffolding adj file for all species from BESST score files (1 pre-scaffolding file for ALL species from ORI !!!)
###
###   INPUT:
###      1- BESST directory
###         (DECOSTAR/DATA_DeCoSTAR/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_ALL/TRIMMOMATIC3/ALL/)
###         (Validation_ADseq/coverage/SCAFFOLDING/BESST-2.2.5/Bowtie2_X_50/50pourc)
###      2- Maximum distance between linked contigs
###         (Ex: 1000000000)
###      3- Links number threshold in scaffolding adjacencies
###         (Ex: 3)
###      4- OUTPUT file
###         (Validation_ADseq/coverage/ADseq/preScaff_ORI_3)
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/preScaff_Aara_50pourc_3)
###      5- Species name
###         (ALL) => If ALL species have to be processed (ORI)
###         (Anopheles_arabiensis)
###
###   OUTPUT:
###      - scaffolding adj file for all species to create instance for ART-DeCo_seq. 1st step, OUTPUT file need to be refined by:
###           => 05b-scaff_adj_file_final.py
###
###   Name: 06a-create_scaff_adj_prefile_BLASTn.py  Author: Yoann Anselmetti
###   Creation date: 2015/12/02                     Last modification: 2017/03/18
###

from sys import argv, stdout
from re import search
from os import close, listdir, path, makedirs
from datetime import datetime
from collections import namedtuple   #New in version 2.6
import errno

stdout.flush()


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise



def rev_ori(ori):
   if ori=="-":
      return "+"
   elif ori=="+":
      return "-"
   else:
      exit("ERROR, orientation should be \"+\" or \"-\" and not: "+ori)



def order_ADJ(ctg1,ctg2,ori1,ori2):
   out_ctg1=""
   out_ctg2=""
   out_ori1=""
   out_ori2=""
   if ctg1 < ctg2:
      out_ctg1=ctg1
      out_ctg2=ctg2
      out_ori1=ori1
      out_ori2=ori2
   else:
      out_ctg1=ctg2
      out_ctg2=ctg1
      out_ori1=rev_ori(ori2)
      out_ori2=rev_ori(ori1)
   return out_ctg1,out_ctg2,out_ori1,out_ori2



def add_ADJ(adj,edge,dict_spe_edge_scaff,score1,score2,links_nb):
   # If the current edge is already in dict_spe_edge_scaff.
   if adj in dict_spe_edge_scaff[species]:
      if dict_spe_edge_scaff[species][adj].vscore<score1:
         dict_spe_edge_scaff[species][adj]=edge
      elif dict_spe_edge_scaff[species][adj].vscore==score1:
         # If vscore are equals take the adj with the best dscore
         if dict_spe_edge_scaff[species][adj].dscore<score2:
            dict_spe_edge_scaff[species][adj]=edge
         elif dict_spe_edge_scaff[species][adj].dscore==score2:
            # If dscore are equals take the adaj with the higher number of links
            if dict_spe_edge_scaff[species][adj].link<int(links_nb):
               dict_spe_edge_scaff[species][adj]=edge
   else:
      dict_spe_edge_scaff[species][adj]=edge
   return dict_spe_edge_scaff



def store_ADJ(BESST_dir,species,dict_spe_edge_scaff):
   for SRX in sorted(listdir(BESST_dir+"/"+species)):
      if SRX!="BESST_output" and SRX!="newADJ_BESST_complete":
         for SRR in sorted(listdir(BESST_dir+"/"+species+"/"+SRX)):
            ctg_scaff_graph_file=""
            # If BAM merged
            if SRR=="BESST_output":
               ctg_scaff_graph_file=BESST_dir+"/"+species+"/"+SRX+"/BESST_output/score_file_pass_1.tsv"
            # If BAM NOT merged
            else:
               ctg_scaff_graph_file=BESST_dir+"/"+species+"/"+SRX+"/"+SRR+"/BESST_output/score_file_pass_1.tsv"

            # print ctg_scaff_graph_file
            scaff_graph=open(ctg_scaff_graph_file,'r')
            for line in scaff_graph:
               r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\n\t]*)\n$',line)
               if r:
                  ctg1=r.group(1)
                  ori1=r.group(2)
                  ctg2=r.group(3)
                  ori2=r.group(4)
                  gap=r.group(5)
                  var_score=r.group(6)
                  disp_score=r.group(7)
                  links_nb=r.group(8)

                  if ctg1!="scf1/ctg1":
                     if int(links_nb)>links_min and float(gap)<=gap_max:
                        score1='{0:.12f}'.format(float(var_score))
                        score2='{0:.12f}'.format(float(disp_score))

                        # ORDER CTG IN ADJ TO HAVE THE SAME ORDER IN PAIRS OF CTG 
                        ctg1,ctg2,ori1,ori2 = order_ADJ(ctg1,ctg2,ori1,ori2)

                        adj=ADJ(species,ctg1,ctg2,ori1,ori2)
                        edge=EDGE(species,ctg1,ctg2,ori1,ori2,float(gap),score1,score2,int(links_nb))
                        # Add adj in dict_spe_edge_scaff if not present or if better score
                        dict_spe_edge_scaff = add_ADJ(adj,edge,dict_spe_edge_scaff,score1,score2,links_nb)

               else:
                  exit("ERROR, the line:\n\t"+line+"\nis incorrectly written in file "+ctg_scaff_graph_file)
            scaff_graph.close()
   return dict_spe_edge_scaff



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   BESST_dir=argv[1]
   gap_max=float(argv[2])
   links_min=int(argv[3])
   OUTPUT_file=argv[4]
   TAG=argv[5]

   OUTPUT_DIR=path.dirname(OUTPUT_file)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR=path.normpath(OUTPUT_DIR)

   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)


   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   EDGE=namedtuple("EDGE",["spe","ctg1","ctg2","ori1","ori2","gap","vscore","dscore","link"])
   ADJ=namedtuple("ADJ",["spe","ctg1","ctg2","ori1","ori2"])


##########
### BROWSE BESST DIRECTORY TO GET SCAFFOLDING ADJACENCIES PROPOSED BY BESST AND STORE IT IN dict_spe_edge_scaff
##########
   dict_spe_edge_scaff={}
   if TAG=="ALL":
      for species in sorted(listdir(BESST_dir)):
         print species
         dict_spe_edge_scaff[species]=dict()
         dict_spe_edge_scaff = store_ADJ(BESST_dir,species,dict_spe_edge_scaff)
   else:
      species=TAG
      print species
      dict_spe_edge_scaff[species]=dict()
      dict_spe_edge_scaff = store_ADJ(BESST_dir,species,dict_spe_edge_scaff)

         

##########
### WRITE OUTPUT_FILE => PRE SCAFFOLDING GRAPH FILE
##########
   output=open(OUTPUT_file,"w")
   output.write("species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tctg1-ctg2_dist\tvscore\tdscore\t#links\n")
   for spe in sorted(dict_spe_edge_scaff):
      for adj in sorted(dict_spe_edge_scaff[spe]):
         output.write(dict_spe_edge_scaff[spe][adj].spe+"\t"+dict_spe_edge_scaff[spe][adj].ctg1+"\t"+dict_spe_edge_scaff[spe][adj].ctg2+"\t"+dict_spe_edge_scaff[spe][adj].ori1+"\t"+dict_spe_edge_scaff[spe][adj].ori2+"\t"+str(dict_spe_edge_scaff[spe][adj].gap)+"\t"+str(dict_spe_edge_scaff[spe][adj].vscore)+"\t"+str(dict_spe_edge_scaff[spe][adj].dscore)+"\t"+str(dict_spe_edge_scaff[spe][adj].link)+"\n")
   output.close()


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))

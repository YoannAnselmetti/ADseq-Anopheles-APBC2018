#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###       Create scaffolding adj file for all species from BESST score files
###
###   INPUT:
###      1- PRE-Scaffolding adjacencies file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/preScaff_Aara_50pourc_3)
###      2- File of association between NEW and OLD CTG/SCAFF minia ID
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/assoCTG_Aara_50pourc)
###      3- OUTPUT PRE-scaffolding adjacencies file modified
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/preScaff_Aara_50pourc_3_modif)
###
###   OUTPUT:
###      - scaffolding adj file for all species to create instance for ART-DeCo_seq. 1st step, OUTPUT file need to be refined by:
###           => 05b-scaff_adj_file_final.py
###
###   Name: 06b-modif_scaff_adj_prefile_BLASTn.py  Author: Yoann Anselmetti
###   Creation date: 2015/12/02                    Last modification: 2017/03/18
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



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   preScaff_file=open(argv[1],'r')
   assoCTG_file=open(argv[2],'r')
   OUTPUT_file=argv[3]

   OUTPUT_DIR=path.dirname(OUTPUT_file)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR=path.normpath(OUTPUT_DIR)

   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   CTG=namedtuple("CTG",["ori","Id"])
   ADJ_ORI=namedtuple("ADJ",["ctg1","ctg2","o1","o2"])
   EDGE=namedtuple("EDGE",["spe","ctg1","ctg2","ori1","ori2","gap","vscore","dscore","link"])


###############################################################################################################
### BROWSE BESST DIRECTORY TO GET SCAFFOLDING ADJACENCIES PROPOSED BY BESST AND STORE IT IN dict_edge_scaff ###
###############################################################################################################
   dict_old_new=dict()
   for line in assoCTG_file:
      r=search('^([^\t]*)\t([^\t\n]*)\n$',line)
      if r:
         new_CTG=r.group(1)
         old_CTG_list=r.group(2)
         if new_CTG!="NEW_CTG_ID":
            # print new_CTG+":"
            for old_CTG in old_CTG_list.split(":"):
               r_ctg=search('^([+-])(.*)',old_CTG)
               if r_ctg:
                  ori=r_ctg.group(1)
                  ctg=r_ctg.group(2)
                  # print "\t"+ctg+" ("+ori+")"
                  dict_old_new[ctg]=CTG(ori,new_CTG)
               else:
                  exit("ERROR in CTG "+old_CTG)
   assoCTG_file.close()

   dict_spe_ADJ=dict()
   output=open(OUTPUT_file,'w')
   for line in preScaff_file:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$',line)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         ori_ctg1=r.group(4)
         ori_ctg2=r.group(5)
         gap=r.group(6)
         vscore=r.group(7)
         dscore=r.group(8)
         links=r.group(9)

         if species=="species":
            output.write(line)
         else:
            if not species in dict_spe_ADJ:
               dict_spe_ADJ[species]=dict()

            if (ctg1 in dict_old_new and ctg2 in dict_old_new):
               new_ctg1=dict_old_new[ctg1]
               new_ctg2=dict_old_new[ctg2]
               if (new_ctg1.Id != new_ctg2.Id):
                  ori1=""
                  if (new_ctg1.ori=="+"):
                     ori1=ori_ctg1
                  elif (new_ctg1.ori=="-"):
                     ori1 = rev_ori(ori_ctg1)
                  else:
                     exit("ERROR, orientation of new_ctg1 "+new_ctg1.Id+" is unknown "+ori_ctg1.ori)
                  ori2=""
                  if (new_ctg2.ori=="+"):
                     ori2=ori_ctg2
                  elif (new_ctg2.ori=="-"):
                     ori2 = rev_ori(ori_ctg2)
                  else:
                     exit("ERROR, orientation of new_ctg2 "+new_ctg2.Id+" is unknown "+ori_ctg2.ori)


                  ctg1,ctg2,ori1,ori2 = order_ADJ(new_ctg1.Id,new_ctg2.Id,ori1,ori2)

                  score1='{0:.12f}'.format(float(vscore))
                  score2='{0:.12f}'.format(float(dscore))
                  adj=ADJ_ORI(ctg1,ctg2,ori1,ori2)
                  edge=EDGE(species,ctg1,ctg2,ori1,ori2,gap,score1,score2,links)

                  if adj in dict_spe_ADJ[species]:
                     if dict_spe_ADJ[species][adj].vscore<score1:
                        dict_spe_ADJ[species][adj]=edge
                     elif dict_spe_ADJ[species][adj].vscore==score1:
                        # If vscore are equals take the adj with the best dscore
                        if dict_spe_ADJ[species][adj].dscore<score2:
                           dict_spe_ADJ[species][adj]=edge
                        elif dict_spe_ADJ[species][adj].dscore==score2:
                           # If dscore are equals take the adj with the higher number of links
                           if dict_spe_ADJ[species][adj].link<int(links):
                              dict_spe_ADJ[species][adj]=edge
                  else:
                     dict_spe_ADJ[species][adj]=edge

   for spe in sorted(dict_spe_ADJ):
      for adj in sorted(dict_spe_ADJ[spe]):
         edge=dict_spe_ADJ[spe][adj]
         output.write(edge.spe+"\t"+edge.ctg1+"\t"+edge.ctg2+"\t"+edge.ori1+"\t"+edge.ori2+"\t"+edge.gap+"\t"+str(edge.vscore)+"\t"+str(edge.dscore)+"\t"+edge.link+"\n")
   output.close()

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))

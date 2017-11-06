#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Get new ADJ predicted by BESST in a complete run not only with score support
###
###   INPUT:
###      1- AGP file BESST results
###         (data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_k50/ERNE-filter/blastn/50pourc/Anopheles_albimanus/BESST_output/pass2/info-pass2.agp)
###      2- CTG file  
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_albimanus/50pourc/CTG_Aalb_50pourc_filt)
###      3- Output file with ADJ predicted by BESST scaffolding algorithm (full run)
###         (data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_k50/ERNE-filter/blastn/50pourc/Anopheles_albimanus/newADJ_BESST_complete)
###      4- Species name
###         (Anopheles_albimanus)
###
###   OUTPUT:
###      - Give adjacencies that are in adjacencies file proposed by DeCo* (param 1) and are inconsistent with adajcencies in reference genome assembly (param 3)
###
###   Name: transform_AGPfile_in_ADJfile_blastn.py    Author: Yoann Anselmetti
###   Creation date: 2017/01/05                      Last modification: 2017/03/20
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


def store_CTGinfo(CTG_file):
   contig_file=open(CTG_file,'r')
   dict_ID_ctg={}
   # When get a contigs pairs (edge scaffolding link) => Get genes that are linked by scaffolding graph with the distance
   for line in contig_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         spe=r.group(1)
         contig=r.group(2)
         contig_size=r.group(3)
         contig_geneNb=r.group(4)
         GF1=r.group(5)
         g1=r.group(6)
         oriG1=r.group(7)
         start_g1=r.group(8)
         GF2=r.group(9)
         g2=r.group(10)
         oriG2=r.group(11)
         stop_g2=r.group(12)

         if spe==select_spe:
            if contig_size=="?":
               print "\n\t=> Contig "+contig+" is not present in FASTA file assembly of species "+spe
            else:
               ctg=CTG(spe,contig,int(contig_size),GF1,g1,oriG1,int(start_g1),GF2,g2,oriG2,int(stop_g2))
               dict_ID_ctg[contig]=ctg
      else:
         exit("ERROR in line "+line+" of file "+contigEXT_file+" !!!")
   contig_file.close()

   return dict_ID_ctg



def parse_AGPfile(agp_file):
   dict_scaff_ctg=dict()
   # Browse edge from scaffolding graph file and store them in "dict_edge_scaff" (1 scaffolding graph file for all species)
   for line in agp_file:
      r1=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r1:
         # print line
         scaff=r1.group(1)
         scaff_start=r1.group(2)
         scaff_end=r1.group(3)
         ID=r1.group(4)
         compo_type=r1.group(5)
         ctg=r1.group(6)
         ctg_start=r1.group(7)
         ctg_end=r1.group(8)
         ctg_ori=r1.group(9)

         if compo_type!="N" and compo_type!="U":
            ori_ctg=CTG_SHORT(ctg_ori,ctg)
            if not scaff in dict_scaff_ctg: 
               dict_scaff_ctg[scaff]=list()
            dict_scaff_ctg[scaff].append(ori_ctg)
   agp_file.close()

   return dict_scaff_ctg



def affect_correct_orientation_ctg1(ori1,ctg1,dict_ID_ctg,CTG_file):
   boolCTG=False
   gene1=""
   oriG1=""
   if ctg1 in dict_ID_ctg:
      boolCTG=True
      if ori1=="+":
         gene1=dict_ID_ctg[ctg1].g2
         oriG1=dict_ID_ctg[ctg1].ori2
      elif ori1=="-":
         gene1=dict_ID_ctg[ctg1].g1
         oriG1=dict_ID_ctg[ctg1].ori1
         if (oriG1=="-"):
            oriG1="+"
         elif (oriG1=="+"):
            oriG1="-"
         else:
            exit("ERROR on gene orientation, should be \"+\" or \"-\" and not: "+oriG1+" !!!")
      else:
         exit("ERROR, orientation of contig, should be \"-\" or \"+\" and not: "+ori1+" !!!")
   else:
      print "\t=> CTG "+ctg1+" has been filtered during filtering of contigs with genes present in gene trees"

   return boolCTG,gene1,oriG1



def affect_correct_orientation_ctg2(ori2,ctg2,dict_ID_ctg,CTG_file):
   boolCTG=False
   gene2=""
   oriG2=""
   if ctg2 in dict_ID_ctg:
      boolCTG=True
      if ori2=="+":
         gene2=dict_ID_ctg[ctg2].g1
         oriG2=dict_ID_ctg[ctg2].ori1
      elif ori2=="-":
         gene2=dict_ID_ctg[ctg2].g2
         oriG2=dict_ID_ctg[ctg2].ori2
         if (oriG2=="-"):
            oriG2="+"
         elif (oriG2=="+"):
            oriG2="-"
         else:
            exit("ERROR on gene orientation, should be \"+\" or \"-\" and not: "+oriG2+" !!!")
      else:
         exit("ERROR, orientation of contig, should be \"-\" or \"+\" and not: "+ori2+" !!!")
   else:
      print "\t=> CTG "+ctg2+" has been filtered during filtering of contigs with genes present in gene trees"

   return boolCTG,gene2,oriG2


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   agp_file=open(argv[1],'r')
   CTG_file=argv[2]
   output_file=open(argv[3],'w')
   select_spe=argv[4]


   output_file.write("#species\tgene1\tgene2\tori1\tori2\n")

   CTG=namedtuple("CTG",["spe","ctg","size","gf1","g1","ori1","start","gf2","g2","ori2","end"])
   CTG_SHORT=namedtuple("CTG_SHORT",["ori","Id"])
   PAIR=namedtuple("PAIR",["ctg1","ctg2","ori1","ori2"])


#########
### STORE CTG INFOS IN dict_ID_ctg
#########
   print "1/ Store CTG infos... ",
   dict_ID_ctg = store_CTGinfo(CTG_file)
   print "DONE"



######################
### Analysed new extant adjacencies predicted by DeCo* in MINIA genome assembly of selected species
   print "2/ Browse contig adjacencies in AGP file to create new adjacencies file predicted by BESST (complete run)... ",
   dict_scaff_ctg = parse_AGPfile(agp_file)
   print "DONE"



   # for scaff in dict_scaff_ctg:
   #    print scaff
   #    for ctg in dict_scaff_ctg[scaff]:
   #       print "\t- "+ctg.Id+" "+ctg.ori



#####################
### Write new adjacencies file predicted by BESST (complete run)
   i=0
   for scaff in dict_scaff_ctg:
      if len(dict_scaff_ctg[scaff])>1:
         # print dict_scaff_ctg[scaff]
         ctg1=""
         ori1=""
         bool1=False
         bool2=False
         for ctg in dict_scaff_ctg[scaff]:
            curr_ctg=ctg.Id
            curr_ori=ctg.ori

            # If curr_ctg is the first CTG of scaff (or at least first with gene)
            if not ctg1:
               if curr_ctg in dict_ID_ctg:
                  ctg1=curr_ctg
                  ori1=curr_ori
               else:
                  print "\t=> CTG "+curr_ctg+" has been filtered during filtering of contigs with genes present in gene trees"

            #If curr_ctg is NOT the first CTG of scaff (is not the first CTG of the scaffold present in dict_ID_ctg)
            else:
               i+=1
               ctg2=curr_ctg
               ori2=curr_ori

               if ctg1!=ctg2:
                  # print ctg
                  # Fill gene1 and oriG1 => If ctg1 is not in dict_ID_ctg then gene1 and oriG1
                  bool1,gene1,oriG1 = affect_correct_orientation_ctg1(ori1,ctg1,dict_ID_ctg,CTG_file)
                  # Fill gene2 and ori2
                  bool2,gene2,oriG2 = affect_correct_orientation_ctg2(ori2,ctg2,dict_ID_ctg,CTG_file)

               else:
                  exit("ERROR, contig "+ctg1+" is adjacent to itself!!!")

               if bool1 and bool2:
                  # Write adjacency given by AGP file in output_file
                  output_file.write(select_spe+"\t"+gene1+"\t"+gene2+"\t"+oriG1+"\t"+oriG2+"\n")
                  ctg1=ctg2
                  ori1=ori2
               else:
                  print "ADJ not considered due to CTG that have been filtered:" 
                  print "\tGENE1: "+ctg1+" "+ori1+" "+gene1+" "+oriG1
                  print "\tGENE2: "+ctg2+" "+ori2+" "+gene2+" "+oriG2


   print "\t=> There are "+str(i)+" new adjacencies."

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
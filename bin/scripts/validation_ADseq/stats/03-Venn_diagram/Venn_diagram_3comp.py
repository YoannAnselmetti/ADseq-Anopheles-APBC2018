#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Compare annotation gene files between original and after mapping of genes on genome assembly computed by minia
###
###   INPUT:
###      1- Venn diagram directory
###         (results/validation_ADseq/spi_20/stats/Venn_diagram/)
###      2- Species tag
###         (Aalb)
###      3- Sequencing tag
###         (50pourc)
###      4- Adjacencies set tag
###         (TP / TFP / FP_TFP / ALL)
###      5- Threshold for ADJ support
###         (0.1 / 0.5 / 0.8)
###      6- boolean indicating if orientation has to be taken into account in adjacency similarity  
###         (T/t/True/true  or F/f/False/false)
###      
###
###   OUTPUT:
###      - output Venn diagram of 
###
###   Name: venn_diagram_3comp.py      Author: Yoann Anselmetti
###   Creation date: 2016/02/15        Last modification: 2017/02/24
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import numpy as np
import errno

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise



def chg_ori(ori):
   if ori=="-":
      return "+"
   elif ori=="+":
      return "-"
   else:
      exit("ERROR, orientation should be \"+\" or \"-\" and not: "+ori)



def clean_ADJ_ori(g1,g2,ori1,ori2,dictADJ,listADJ):
   if g1 < g2:
      adj=ADJ_ORI(g1,g2,ori1,ori2)
      listADJ.append(adj)
      if not g1 in dictADJ:
         dictADJ[g1]=list()
      dictADJ[g1].append(adj)
   else:
      g1_modif=g2
      g2_modif=g1
      ori1_modif=chg_ori(ori2)
      ori2_modif=chg_ori(ori1)
      adj=ADJ_ORI(g1_modif,g2_modif,ori1_modif,ori2_modif)
      listADJ.append(adj)
      if not g1_modif in dictADJ:
         dictADJ[g1_modif]=list()
      dictADJ[g1_modif].append(adj)



def clean_ADJ(g1,g2,dictADJ,listADJ):
   if g1 < g2:
      adj=ADJ(g1,g2)
      listADJ.append(adj)
      if not g1 in dictADJ:
         dictADJ[g1]=list()
      dictADJ[g1].append(adj)
   else:
      g1_modif=g2
      g2_modif=g1
      adj=ADJ(g1_modif,g2_modif)
      listADJ.append(adj)
      if not g1_modif in dictADJ:
         dictADJ[g1_modif]=list()
      dictADJ[g1_modif].append(adj)



def file_dictADJ(file_path,dictADJ,listADJ,bool_ori):
   adj_nb=0
   file=open(file_path,'r')
   for adj in file:
      if bool_ori:
         r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",adj)
         if r:
            g1=r.group(1)
            g2=r.group(2)
            ori1=r.group(3)
            ori2=r.group(4)

            adj_nb+=1
            clean_ADJ_ori(g1,g2,ori1,ori2,dictADJ,listADJ)

         else:
            exit("ERROR, file:\n\t"+file_path+"\nis incorrectly writtten !!!")

      else:
         r=search("^([^\t]*)\t([^\t\n]*)\n$",adj)
         if r:
            g1=r.group(1)
            g2=r.group(2)

            adj_nb+=1
            clean_ADJ(g1,g2,dictADJ,listADJ)

         else:
            exit("ERROR, file:\n\t"+file_path+"\nis incorrectly writtten !!!")

   file.close()
   # print "ADJ: "+str(adj_nb)


def filter_commonADJ(dictA,dictB,dictAB,listA,listB,listAB):
   adj=0
   for gene in dictA:
      if gene in dictB:
         for adj1 in dictA[gene]:
            for adj2 in dictB[gene]:
               if adj1==adj2:
                  if not gene in dictAB:
                     dictAB[gene]=list()
                  dictAB[gene].append(adj1)
                  listAB.append(adj1)
                  if adj1 in listA:
                     listA.remove(adj1)
                  if adj1 in listB:
                     listB.remove(adj1)


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   venn_dir=argv[1]
   spe_tag=argv[2]
   seq_tag=argv[3]
   adj_tag=argv[4]
   support=argv[5]
   orientation=argv[6]
   OUTPUT_DIR=argv[7]


   file1=venn_dir+"/"+spe_tag+"/"+seq_tag+"/"+adj_tag+"/ADseq_LINEAR_"+support
   file2=venn_dir+"/"+spe_tag+"/"+seq_tag+"/"+adj_tag+"/AD_LINEAR_"+support
   file3=venn_dir+"/"+spe_tag+"/"+seq_tag+"/"+adj_tag+"/BESST_complete"


   bool_ori=False
   if orientation=="T" or orientation=="True" or orientation=="true" or orientation=="t" or orientation=="Y" or orientation=="Yes" or orientation=="yes" or orientation=="y":
      bool_ori=True
   elif orientation=="F" or orientation=="False" or orientation=="false" or orientation=="f" or orientation=="N" or orientation=="No" or orientation=="no" or orientation=="n":
      bool_ori=False
   else:
      exit("ERROR, parameter 4 for orientation boolean should be equal to \"true\" or \"false\" and NOT: "+orientation+" !!!")

   ADJ=namedtuple("ADJ",["g1","g2"])
   ADJ_ORI=namedtuple("ADJ",["g1","g2","o1","o2"])



   #########
   ### Store new adjacencies list of INPUT files in dict and list   
   #########
   dict1=dict()
   list1=list()
   print "Store adjacencies of file: \""+file1+"\" in dict1 and list1",
   stdout.flush()
   file_dictADJ(file1,dict1,list1,bool_ori)
   print "DONE"

   dict2=dict()
   list2=list()
   print "Store adjacencies of file: \""+file2+"\" in dict2 and list2",
   stdout.flush()
   file_dictADJ(file2,dict2,list2,bool_ori)
   print "DONE"

   dict3=dict()
   list3=list()
   print "Store adjacencies of file: \""+file3+"\" in dict3 and list3",
   stdout.flush()
   file_dictADJ(file3,dict3,list3,bool_ori)
   print "DONE\n"



   adj1_nb=0
   for g in dict1:
      adj1_nb+=len(dict1[g])
   print "\nlist1 size: "+str(len(list1))+" & dict1 size: "+str(adj1_nb)
   adj2_nb=0
   for g in dict2:
      adj2_nb+=len(dict2[g])
   print "list2 size: "+str(len(list2))+" & dict2 size: "+str(adj2_nb)
   adj3_nb=0
   for g in dict3:
      adj3_nb+=len(dict3[g])
   print "list3 size: "+str(len(list3))+" & dict3 size: "+str(adj3_nb)+"\n"



   # Venn diagramm with three elements
   ###########
   ### Determine common new adajcnecies between the three files to draw Venn diagram.
   ###########
   dict12=dict()
   list12=list()
   print "\nCommon ADj between files: \""+file1+"\" AND \""+file2+"\"",
   stdout.flush()
   filter_commonADJ(dict1,dict2,dict12,list1,list2,list12)
   print "DONE"
   print "\t=> list1 size: "+str(len(list1))
   print "\t=> list2 size: "+str(len(list2))
   print "\t=> list12 size: "+str(len(list12))

   dict13=dict()
   list13=list()
   print "Common ADj between files: \""+file1+"\" AND \""+file3+"\"",
   stdout.flush()
   filter_commonADJ(dict1,dict3,dict13,list1,list3,list13)
   print "DONE"
   print "\t=> list1 size: "+str(len(list1))
   print "\t=> list3 size: "+str(len(list3))
   print "\t=> list13 size: "+str(len(list13))

   dict23=dict()
   list23=list()
   print "Common ADj between files: \""+file2+"\" AND \""+file3+"\"",
   stdout.flush()
   filter_commonADJ(dict2,dict3,dict23,list2,list3,list23)
   print "DONE"
   print "\t=> list2 size: "+str(len(list2))
   print "\t=> list3 size: "+str(len(list3))
   print "\t=> list23 size: "+str(len(list23))


   list123=list()
   listcomm=list12[:]
   for adj in listcomm:
      if adj in list13 and adj in list23:
         list123.append(adj)
         list12.remove(adj)
         list13.remove(adj)
         list23.remove(adj)



   # Print section new adjacencies number
   A=len(list1)
   B=len(list2)
   C=len(list3)

   AB=len(list12)
   AC=len(list13)
   BC=len(list23)

   ABC=len(list123)

   print "Venn diagram components:"
   print "\t- 1: "+str(A)
   print "\t- 2: "+str(B)
   print "\t- 3: "+str(C)
   print "\t- 1 & 2: "+str(AB)   
   print "\t- 1 & 3: "+str(AC)
   print "\t- 2 & 3: "+str(BC)
   print "\t- 1 & 2 & 3: "+str(ABC)


   adj_tot=A+B+C+AB*2+AC*2+BC*2+ABC*3
   print "There are "+str(adj_tot)+" adjacencies in TOTAL in the 3 INPUT files!!!"


   output_dir=OUTPUT_DIR+"/"+spe_tag+"/"+seq_tag+"/"+support
   mkdir_p(output_dir)
   if bool_ori:
      fig_name=output_dir+"/Venn_diagram_"+adj_tag+"_+ori.pdf"
   else:
      fig_name=output_dir+"/Venn_diagram_"+adj_tag+"-ori.pdf"



   v=venn3(subsets=(A,B,AB,C,AC,BC,ABC),set_labels=('ADseq', 'AD', 'BESST'))

   for text in v.set_labels:
      text.set_fontsize(20)

   # for text in v.subset_labels:
   #    text.set_fontsize(18)

   # if bool_ori:
   #    plt.title("Venn diagram: spe: "+spe_tag+", seq: "+seq_tag+", support: "+support+", orientation: True AND adj set: "+adj_tag)
   # else:
   #    plt.title("Venn diagram: spe: "+spe_tag+", seq: "+seq_tag+", support: "+support+", orientation: False AND adj set: "+adj_tag)

   plt.title(adj_tag,size=30)
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
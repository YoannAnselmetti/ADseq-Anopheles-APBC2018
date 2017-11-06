#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Generate plot of Precision and Recall for validation of AD-seq
###
###   INPUT:
###      1- Species name
###         (Anopheles_albimanus / Anopheles_arabiensis / Anopheles_dirus)
###      2- Species TAG
###         (Aalb / Aara / Adir)
###      3- Reads sampling TAG_ID
###         (50pourc / ALL)
###      4- INPUT directory containing all input files
###         (results/validation_ADseq/spi_default/stats)
###         (results/validation_ADseq/spi_20/stats)
###      5- OUTPUT directory where plot will be stored
###         (figures/precision_recall/spi_default)
###         (figures/precision_recall/spi_20)
###
###   OUTPUT files:
###      - Distrbution of Precision and Recall statistics for INPUT species
###
###   Name: plot_Precision_Recall.py         Author: Yoann Anselmetti
###   Creation date: 2017/01/20              Last modification: 2017/02/25
###


from sys import argv, stdout
from re import search
from os import close, path, makedirs, remove

import errno
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise



def plot_recall_precision(list_supp_adseq,list_supp_ad,list_rec_adseq,list_rec_ad,list_prec_adseq,list_prec_ad,rec_besstC,prec_besstC,OUTPUT_DIR):

   fig = plt.figure()
   ax = fig.add_subplot(1,1,1)

   #Plot analytic solution
   ax.plot(list_supp_adseq, list_rec_adseq, color='g', label="ADseq", marker='o',markersize=10)
   ax.plot(list_supp_ad, list_rec_ad, color='g', label="AD", marker='^',markersize=10)
   ax.plot((0,1),(rec_besstC,rec_besstC), color='g', label="BESST(complete)", marker='s',markersize=10)

   #Plot simulation
   ax.plot(list_supp_adseq, list_prec_adseq, color='r', marker='o',markersize=10)
   ax.plot(list_supp_ad, list_prec_ad, color='r', marker='^',markersize=10)
   ax.plot((0,1),(prec_besstC,prec_besstC), color='r', marker='s',markersize=10)
   
   # Artist for the legend
   adseqArtist=plt.Line2D((0,1),(0,0), color='k', marker='o',markersize=10)
   adArtist=plt.Line2D((0,1),(0,0), color='k', marker='^',markersize=10)
   besstCArtist=plt.Line2D((0,1),(0,0), color='k', marker='s',markersize=10)
   recArtist=plt.Line2D((0,1),(0,0), color='g',markersize=10)
   precArtist=plt.Line2D((0,1),(0,0), color='r',markersize=10)

   # Create legend from custom artist/label lists
   ax.legend([adseqArtist,adArtist,besstCArtist,precArtist,recArtist],['ADseq','AD','BESST','Precision','Recall'],numpoints=1,loc="center")

   plt.title(species)
   plt.xlabel("Adjacency support threshold")
   plt.xticks(np.arange(0, 1.1, 0.1))
   plt.xlim([-0.01,1.01])

   plt.ylabel("%")
   plt.yticks(np.arange(0, 101, 10))
   plt.ylim([-1,101])
   
   fig_name=OUTPUT_DIR+"/distrib_rec_prec_"+tag+"_"+reads+".pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()




def read_file_and_store_stats_DeCoSTAR(stats_file,list_prec,list_rec,list_supp):
   statsFile=open(stats_file,"r")
   for line in statsFile:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', line)
      if r:
         filt=r.group(1)
         predicted=r.group(2)
         find=r.group(3)
         FN=r.group(4)
         TP=r.group(5)
         FP=r.group(6)
         CFP=r.group(7)
         CFP_FP=r.group(8)
         recall=r.group(9)
         precision=r.group(10)
         recallCFP=r.group(11)
         precisionCFP=r.group(12)

         if filt!="filter" and not "LIN" in filt and recallCFP!="NA":
            list_supp.append(float(filt))
            list_prec.append(float(precisionCFP))
            list_rec.append(float(recallCFP))
   statsFile.close()



def read_file_and_store_stats_BESST(stats_file,list_stats_besst):
   statsFile=open(stats_file,"r")
   for line in statsFile:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', line)
      if r:
         filt=r.group(1)
         predicted=r.group(2)
         find=r.group(3)
         FN=r.group(4)
         TP=r.group(5)
         FP=r.group(6)
         CFP=r.group(7)
         CFP_FP=r.group(8)
         recall=r.group(9)
         precision=r.group(10)
         recallCFP=r.group(11)
         precisionCFP=r.group(12)


         if filt=="COMPLET" and precisionCFP!="NA" and recallCFP!="NA":
            list_stats_besst.append(float(precisionCFP))
            list_stats_besst.append(float(recallCFP))
   statsFile.close()

# Recovery of input parameters
species=argv[1]
tag=argv[2]
reads=argv[3]
inputDIR=argv[4]
outputDIR=argv[5]


OUTPUT_DIR=outputDIR+"/"+reads

mkdir_p(OUTPUT_DIR)


adseq=inputDIR+"/"+tag+"_"+reads+"_ADseq"
ad=inputDIR+"/"+tag+"_"+reads+"_AD"
besst=inputDIR+"/"+tag+"_"+reads+"_BESST"
# ragout=inputDIR+"/"+tag+"_"+reads+"_BESST"



list_supp_adseq=list()
list_supp_ad=list()

list_prec_adseq=list()
list_rec_adseq=list()

list_prec_ad=list()
list_rec_ad=list()


list_stats_besst=list()

read_file_and_store_stats_DeCoSTAR(adseq,list_prec_adseq,list_rec_adseq,list_supp_adseq)
read_file_and_store_stats_DeCoSTAR(ad,list_prec_ad,list_rec_ad,list_supp_ad)
read_file_and_store_stats_BESST(besst,list_stats_besst)

prec_besstC=0.0
rec_besstC=0.0
if list_stats_besst:
   prec_besstC=list_stats_besst[0]
   rec_besstC=list_stats_besst[1]


plot_recall_precision(list_supp_adseq,list_supp_ad,list_rec_adseq,list_rec_ad,list_prec_adseq,list_prec_ad,rec_besstC,prec_besstC,OUTPUT_DIR)
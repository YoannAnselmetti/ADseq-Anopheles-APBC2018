#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Generate plot of Precision in function of Recall for validation of AD-seq with others scaffolding methods
###
###   INPUT:
###      1- Reads sampling TAG_ID
###         (50pourc / ALL)
###      2- INPUT directory containing all input files
###         (results/validation_ADseq/spi_default/stats)
###         (results/validation_ADseq/spi_20/stats)
###      3- OUTPUT directory where plot will be stored
###         (figures/precision_recall/spi_default)
###         (figures/precision_recall/spi_20)
###
###   Name: plot_Precision_ON_Recall.py      Author: Yoann Anselmetti
###   Creation date: 2017/01/22              Last modification: 2017/01/27
###


from sys import argv, stdout
from re import search
from os import close, path, makedirs, remove
import errno
import fileinput
import subprocess
import itertools
import matplotlib
matplotlib.use('Agg')
 
import numpy as np 
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar

stdout.flush()

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise



def plot_recall_on_precision(rec_adseq_01_Aalb,prec_adseq_01_Aalb,rec_adseq_05_Aalb,prec_adseq_05_Aalb,rec_adseq_08_Aalb,prec_adseq_08_Aalb,rec_ad_01_Aalb,prec_ad_01_Aalb,rec_ad_05_Aalb,prec_ad_05_Aalb,rec_ad_08_Aalb,prec_ad_08_Aalb,rec_besstC_Aalb,prec_besstC_Aalb,rec_adseq_01_Aara,prec_adseq_01_Aara,rec_adseq_05_Aara,prec_adseq_05_Aara,rec_adseq_08_Aara,prec_adseq_08_Aara,rec_ad_01_Aara,prec_ad_01_Aara,rec_ad_05_Aara,prec_ad_05_Aara,rec_ad_08_Aara,prec_ad_08_Aara,rec_besstC_Aara,prec_besstC_Aara,rec_adseq_01_Adir,prec_adseq_01_Adir,rec_adseq_05_Adir,prec_adseq_05_Adir,rec_adseq_08_Adir,prec_adseq_08_Adir,rec_ad_01_Adir,prec_ad_01_Adir,rec_ad_05_Adir,prec_ad_05_Adir,rec_ad_08_Adir,prec_ad_08_Adir,rec_besstC_Adir,prec_besstC_Adir,OUTPUT_DIR):

   fig = plt.figure()
   ax = fig.add_subplot(1,1,1)


   # set_markersize(1.5)

   ### FOR SPECIES ALBIMANUS
   if rec_adseq_01_Aalb>0.0 and prec_adseq_01_Aalb>0.0:
      ax.scatter(rec_adseq_01_Aalb, prec_adseq_01_Aalb, color=(0.5,0,1), label="ADseq(>=0.1)", marker='o',s=100)
   if rec_adseq_05_Aalb>0.0 and prec_adseq_05_Aalb>0.0:
      ax.scatter(rec_adseq_05_Aalb, prec_adseq_05_Aalb, color=(0.75,0,1), label="ADseq(>=0.5)", marker='o',s=100)
   if rec_adseq_05_Aalb>0.0 and prec_adseq_05_Aalb>0.0:
      ax.scatter(rec_adseq_08_Aalb, prec_adseq_08_Aalb, color=(1,0,1), label="ADseq(>=0.8)", marker='o',s=100)

   if rec_ad_01_Aalb>0.0 and prec_ad_01_Aalb>0.0:
      ax.scatter(rec_ad_01_Aalb, prec_ad_01_Aalb, color=(0.5,0,0), label="AD(>=0.1)", marker='o',s=100)
   if rec_ad_05_Aalb>0.0 and prec_ad_05_Aalb>0.0:
      ax.scatter(rec_ad_05_Aalb, prec_ad_05_Aalb, color=(0.75,0,0), label="AD(>=0.5)", marker='o',s=100)
   if rec_ad_08_Aalb>0.0 and prec_ad_08_Aalb>0.0:
      ax.scatter(rec_ad_08_Aalb, prec_ad_08_Aalb, color=(1,0,0), label="AD(>=0.8)", marker='o',s=100)

   if rec_besstC_Aalb>0.0 and prec_besstC_Aalb>0.0:
      ax.scatter(rec_besstC_Aalb, prec_besstC_Aalb, color='b', label="BESST", marker='o',s=100)


   ### FOR SPECIES ARABIENSIS
   if rec_adseq_01_Aara>0.0 and prec_adseq_01_Aara>0.0:
      ax.scatter(rec_adseq_01_Aara, prec_adseq_01_Aara, color=(0.5,0,1), label="ADseq(>=0.1)", marker='s',s=100)
   if rec_adseq_05_Aara>0.0 and prec_adseq_05_Aara>0.0:
      ax.scatter(rec_adseq_05_Aara, prec_adseq_05_Aara, color=(0.75,0,1), label="ADseq(>=0.5)", marker='s',s=100)
   if rec_adseq_05_Aara>0.0 and prec_adseq_05_Aara>0.0:
      ax.scatter(rec_adseq_08_Aara, prec_adseq_08_Aara, color=(1,0,1), label="ADseq(>=0.8)", marker='s',s=100)

   if rec_ad_01_Aara>0.0 and prec_ad_01_Aara>0.0:
      ax.scatter(rec_ad_01_Aara, prec_ad_01_Aara, color=(0.5,0,0), label="AD(>=0.1)", marker='s',s=100)
   if rec_ad_05_Aara>0.0 and prec_ad_05_Aara>0.0:
      ax.scatter(rec_ad_05_Aara, prec_ad_05_Aara, color=(0.75,0,0), label="AD(>=0.5)", marker='s',s=100)
   if rec_ad_08_Aara>0.0 and prec_ad_08_Aara>0.0:
      ax.scatter(rec_ad_08_Aara, prec_ad_08_Aara, color=(1,0,0), label="AD(>=0.8)", marker='s',s=100)

   if rec_besstC_Aara>0.0 and prec_besstC_Aara>0.0:
      ax.scatter(rec_besstC_Aara, prec_besstC_Aara, color='b', label="BESST", marker='s',s=100)


   ### FOR SPECIES DIRUS
   if rec_adseq_01_Adir>0.0 and prec_adseq_01_Adir>0.0:
      ax.scatter(rec_adseq_01_Adir, prec_adseq_01_Adir, color=(0.5,0,1), label="ADseq(>=0.1)", marker='^',s=100)
   if rec_adseq_05_Adir>0.0 and prec_adseq_05_Adir>0.0:
      ax.scatter(rec_adseq_05_Adir, prec_adseq_05_Adir, color=(0.75,0,1), label="ADseq(>=0.5)", marker='^',s=100)
   if rec_adseq_05_Adir>0.0 and prec_adseq_05_Adir>0.0:
      ax.scatter(rec_adseq_08_Adir, prec_adseq_08_Adir, color=(1,0,1), label="ADseq(>=0.8)", marker='^',s=100)

   if rec_ad_01_Adir>0.0 and prec_ad_01_Adir>0.0:
      ax.scatter(rec_ad_01_Adir, prec_ad_01_Adir, color=(0.5,0,0), label="AD(>=0.1)", marker='^',s=100)
   if rec_ad_05_Adir>0.0 and prec_ad_05_Adir>0.0:
      ax.scatter(rec_ad_05_Adir, prec_ad_05_Adir, color=(0.75,0,0), label="AD(>=0.5)", marker='^',s=100)
   if rec_ad_08_Adir>0.0 and prec_ad_08_Adir>0.0:
      ax.scatter(rec_ad_08_Adir, prec_ad_08_Adir, color=(1,0,0), label="AD(>=0.8)", marker='^',s=100)

   if rec_besstC_Adir>0.0 and prec_besstC_Adir>0.0:
      ax.scatter(rec_besstC_Adir, prec_besstC_Adir, color='b', label="BESST", marker='^',s=100)




   ##########
   ### IMSHOW are not available in legend !!!
   ##########
   # m = np.zeros((1,3))
   # for i in range(3):
   #    m[0,i] = i
   # reds = [(0.5,0,0), (0.75,0,0), (1,0,0)]
   # greens = [(0,0.5,0), (0,0.75,0), (0,1,0)]
   # cm_adseq = matplotlib.colors.LinearSegmentedColormap.from_list('adseq', reds)
   # cm_ad = matplotlib.colors.LinearSegmentedColormap.from_list('adseq', greens)

   # fig_adseq, ax_adseq=plt.subplots(figsize=(2,1))
   # fig_ad, ax_ad=plt.subplots(figsize=(2,1))

   # adseqArtist= ax_adseq.imshow(m, cmap=cm_adseq, aspect=2), plt.yticks(np.arange(0)), plt.xticks(np.arange(-0.5,3.5,1.5), [0.1,0.5,0.8])
   # adArtist= ax_ad.imshow(m, cmap=cm_ad, aspect=2), plt.yticks(np.arange(0)), plt.xticks(np.arange(-0.5,3.5,1.5), [0.1,0.5,0.8])



   # Artist for the legend
   AalbArtist=plt.Line2D((0,1),(0,0), color='k', marker='o')
   AaraArtist=plt.Line2D((0,1),(0,0), color='k', marker='s')
   AdirArtist=plt.Line2D((0,1),(0,0), color='k', marker='^')


   adseq_patch = mpatches.Patch(color=(0.75,0,0))
   ad_patch = mpatches.Patch(color=(0,0.75,0))
   besst_patch = mpatches.Patch(color='b')


   # ax.legend(scatterpoints=1,loc="lower right")
   ax.legend([adseq_patch,ad_patch,besst_patch,AalbArtist,AaraArtist,AdirArtist],['ADseq','AD','BESST','A.alb','A.ara','A.dir'],numpoints=1,loc="lower right")


   # plt.title(species)
   plt.xlabel("Recall (%)")
   plt.xticks(np.arange(Xmin, Xmax+1.0, Xscale))
   plt.xlim([Xmin-0.1,Xmax+0.1])

   plt.ylabel("Precision (%)")
   plt.yticks(np.arange(Ymin, Ymax+1.0, Yscale))
   plt.ylim([Ymin-0.1,Ymax+0.1])
   # plt.yticks(np.arange(75, 101, 5))
   # plt.ylim([75,100.1])

   fig_name=OUTPUT_DIR+"/distrib_prec_on_rec_"+reads+".pdf"
   
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()



# ####
# ## PART TO DRAW Gradient color legend (to add in the figure)
# ####
# fig = plt.figure()
# m = np.zeros((1,3))
# for i in range(3):
#     m[0,i] = i
# print m
# reds = [(0.5,0,0), (0.75,0,0), (1,0,0)]
# greens = [(0,0.5,0), (0,0.75,0), (0,1,0)]
# purples = [(0.5,0,1), (0.75,0,1), (1,0,1)]
# cm_R = matplotlib.colors.LinearSegmentedColormap.from_list('R', reds)
# cm_G = matplotlib.colors.LinearSegmentedColormap.from_list('G', greens)
# cm_P = matplotlib.colors.LinearSegmentedColormap.from_list('P', purples)

# plt.imshow(m, cmap=cm_P, aspect=2)
# plt.yticks(np.arange(0))
# plt.xticks(np.arange(-0.5,3.5,1.5), [0.1,0.5,0.8])

# plt.savefig("gradient_color_legend.pdf",format="pdf")
# exit()





def read_file_and_store_stats_DeCoSTAR(stats_file,list_prec,list_rec):
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

         if "LIN" in filt and recallCFP!="NA":
            if "0.1" in filt:
               list_rec.append(float(recallCFP))
               list_prec.append(float(precisionCFP))
            if "0.5" in filt:
               list_rec.append(float(recallCFP))
               list_prec.append(float(precisionCFP))
            if "0.8" in filt:
               list_rec.append(float(recallCFP))
               list_prec.append(float(precisionCFP))
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
            list_stats_besst.append(float(recallCFP))
            list_stats_besst.append(float(precisionCFP))
   statsFile.close()



# Recovery of input parameters
reads=argv[1]
inputDIR=argv[2]
outputDIR=argv[3]

Xmin=0
Xmax=100
Xscale=10
Ymin=50
Ymax=100
Yscale=10

OUTPUT_DIR=outputDIR+"/"+reads


adseq_Aalb=inputDIR+"/Aalb_"+reads+"_ADseq"
ad_Aalb=inputDIR+"/Aalb_"+reads+"_AD"
besst_Aalb=inputDIR+"/Aalb_"+reads+"_BESST"

adseq_Aara=inputDIR+"/Aara_"+reads+"_ADseq"
ad_Aara=inputDIR+"/Aara_"+reads+"_AD"
besst_Aara=inputDIR+"/Aara_"+reads+"_BESST"

adseq_Adir=inputDIR+"/Adir_"+reads+"_ADseq"
ad_Adir=inputDIR+"/Adir_"+reads+"_AD"
besst_Adir=inputDIR+"/Adir_"+reads+"_BESST"



############
### FOR ALBIMANUS
############
list_prec_adseq_Aalb=list()
list_rec_adseq_Aalb=list()
list_prec_ad_Aalb=list()
list_rec_ad_Aalb=list()
list_besst_Aalb=list()
read_file_and_store_stats_DeCoSTAR(adseq_Aalb,list_prec_adseq_Aalb,list_rec_adseq_Aalb)
read_file_and_store_stats_DeCoSTAR(ad_Aalb,list_prec_ad_Aalb,list_rec_ad_Aalb)
read_file_and_store_stats_BESST(besst_Aalb,list_besst_Aalb)


rec_adseq_01_Aalb=list_rec_adseq_Aalb[0]
rec_adseq_05_Aalb=list_rec_adseq_Aalb[1]
rec_adseq_08_Aalb=list_rec_adseq_Aalb[2]
prec_adseq_01_Aalb=list_prec_adseq_Aalb[0]
prec_adseq_05_Aalb=list_prec_adseq_Aalb[1]
prec_adseq_08_Aalb=list_prec_adseq_Aalb[2]

rec_ad_01_Aalb=list_rec_ad_Aalb[0]
rec_ad_05_Aalb=list_rec_ad_Aalb[1]
rec_ad_08_Aalb=list_rec_ad_Aalb[2]
prec_ad_01_Aalb=list_prec_ad_Aalb[0]
prec_ad_05_Aalb=list_prec_ad_Aalb[1]
prec_ad_08_Aalb=list_prec_ad_Aalb[2]

rec_besstC_Aalb=list_besst_Aalb[0]
prec_besstC_Aalb=list_besst_Aalb[1]



############
### FOR ARABIENSIS
############
list_prec_adseq_Aara=list()
list_rec_adseq_Aara=list()
list_prec_ad_Aara=list()
list_rec_ad_Aara=list()
list_besst_Aara=list()
read_file_and_store_stats_DeCoSTAR(adseq_Aara,list_prec_adseq_Aara,list_rec_adseq_Aara)
read_file_and_store_stats_DeCoSTAR(ad_Aara,list_prec_ad_Aara,list_rec_ad_Aara)
read_file_and_store_stats_BESST(besst_Aara,list_besst_Aara)


rec_adseq_01_Aara=list_rec_adseq_Aara[0]
rec_adseq_05_Aara=list_rec_adseq_Aara[1]
rec_adseq_08_Aara=list_rec_adseq_Aara[2]
prec_adseq_01_Aara=list_prec_adseq_Aara[0]
prec_adseq_05_Aara=list_prec_adseq_Aara[1]
prec_adseq_08_Aara=list_prec_adseq_Aara[2]

rec_ad_01_Aara=list_rec_ad_Aara[0]
rec_ad_05_Aara=list_rec_ad_Aara[1]
rec_ad_08_Aara=list_rec_ad_Aara[2]
prec_ad_01_Aara=list_prec_ad_Aara[0]
prec_ad_05_Aara=list_prec_ad_Aara[1]
prec_ad_08_Aara=list_prec_ad_Aara[2]

rec_besstC_Aara=0.0
prec_besstC_Aara=0.0
if list_besst_Aara:
   rec_besstC_Aara=list_besst_Aara[0]
   prec_besstC_Aara=list_besst_Aara[1]

############
### FOR DIRUS
############
list_prec_adseq_Adir=list()
list_rec_adseq_Adir=list()
list_prec_ad_Adir=list()
list_rec_ad_Adir=list()
list_besst_Adir=list()
read_file_and_store_stats_DeCoSTAR(adseq_Adir,list_prec_adseq_Adir,list_rec_adseq_Adir)
read_file_and_store_stats_DeCoSTAR(ad_Adir,list_prec_ad_Adir,list_rec_ad_Adir)
read_file_and_store_stats_BESST(besst_Adir,list_besst_Adir)


rec_adseq_01_Adir=list_rec_adseq_Adir[0]
rec_adseq_05_Adir=list_rec_adseq_Adir[1]
rec_adseq_08_Adir=list_rec_adseq_Adir[2]
prec_adseq_01_Adir=list_prec_adseq_Adir[0]
prec_adseq_05_Adir=list_prec_adseq_Adir[1]
prec_adseq_08_Adir=list_prec_adseq_Adir[2]

rec_ad_01_Adir=list_rec_ad_Adir[0]
rec_ad_05_Adir=list_rec_ad_Adir[1]
rec_ad_08_Adir=list_rec_ad_Adir[2]
prec_ad_01_Adir=list_prec_ad_Adir[0]
prec_ad_05_Adir=list_prec_ad_Adir[1]
prec_ad_08_Adir=list_prec_ad_Adir[2]

rec_besstC_Adir=list_besst_Adir[0]
prec_besstC_Adir=list_besst_Adir[1]



plot_recall_on_precision(rec_adseq_01_Aalb,prec_adseq_01_Aalb,rec_adseq_05_Aalb,prec_adseq_05_Aalb,rec_adseq_08_Aalb,prec_adseq_08_Aalb,rec_ad_01_Aalb,prec_ad_01_Aalb,rec_ad_05_Aalb,prec_ad_05_Aalb,rec_ad_08_Aalb,prec_ad_08_Aalb,rec_besstC_Aalb,prec_besstC_Aalb,rec_adseq_01_Aara,prec_adseq_01_Aara,rec_adseq_05_Aara,prec_adseq_05_Aara,rec_adseq_08_Aara,prec_adseq_08_Aara,rec_ad_01_Aara,prec_ad_01_Aara,rec_ad_05_Aara,prec_ad_05_Aara,rec_ad_08_Aara,prec_ad_08_Aara,rec_besstC_Aara,prec_besstC_Aara,rec_adseq_01_Adir,prec_adseq_01_Adir,rec_adseq_05_Adir,prec_adseq_05_Adir,rec_adseq_08_Adir,prec_adseq_08_Adir,rec_ad_01_Adir,prec_ad_01_Adir,rec_ad_05_Adir,prec_ad_05_Adir,rec_ad_08_Adir,prec_ad_08_Adir,rec_besstC_Adir,prec_besstC_Adir,OUTPUT_DIR)
#!/usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Read saffolding adjacencies file output by BESST and plot distribution of scores
###
###   INPUT:
###      1 - File with list of values to plot histogram
###         (data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3)
###         (data/data_DeCoSTAR/scaff_BESST_DeCoSTAR)
###      2 - OUTPUT path where graph will be writen
###         (figures/besst_score)
###      3 - Text for the legend of the graph (To write between "")
###         (RAW)
###         (DeCoSTAR)
###
###   OUTPUT:
###      - Distribution graph of BESST scores
###
###   Name: plot_hist_scaff_BESST.py         Author: Yoann Anselmetti
###   Creation date: 2017/02/02              Last modification: 2017/10/24
###


from sys import argv, exit
from re import search
from os import close, path, makedirs, listdir, mkdir
from datetime import datetime

import errno
import operator

import itertools
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import colors






def hist_all(list1,list2,list3,OUTPUT,tag):

   colors=['red','blue','purple']
   labels=["link dispersity score","link variation score","mean link score"]
   hist_param=dict(color=colors,label=labels,bins=10,range=(0,1))

   plt.hist((list1,list2,list3),**hist_param)

   # plt.title("Distribution of BESST adjacency scores")

   plt.legend(prop={'size':12})

   plt.xlabel("Score value")
   plt.ylabel("Frequency")

   plt.xlim(xmin=0, xmax=1)
   plt.ylim(ymin=0)

   plt.xticks(np.arange(0, 1.1, .1))

   fig_name=OUTPUT+"/hist_BESST_scores_ALL_"+tag+".pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()





def hist_spe(dico,OUTPUT,tag):
   # colors=['red','blue','yellow','purple','orange','green','brown','black','grey','coral','fuchsia','crimson','peru','olive','darkviolet','cyan']
   colors=['red','blue','yellow','purple','orange','green','brown','black','grey','coral','fuchsia','crimson','peru','olive']
   labels=[key for key, value in dico.iteritems()]
   hist_param=dict(color=colors,label=labels,bins=10,range=(0,1))

   print labels
   # (value for key, value in dico.iteritems())

   plt.hist([value for key, value in dico.iteritems()],**hist_param)

   # plt.title(legend)

   # plt.legend()

   plt.xlabel("Score value interval")
   plt.ylabel("Frequency")

   plt.xlim(xmin=0, xmax=1)
   plt.ylim(ymin=0)

   plt.xticks(np.arange(0, 1.1, .1))

   fig_name=OUTPUT+"/hist_BESST_scores_SPE_"+tag+".pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()






def hist_spe_norm(dico,OUTPUT,tag):

   dico_refine=dict()
   dict_spe_tag=dict()
   for spe in dico:
      list_score=dico[spe]
      a=np.asarray(list_score)

      hist,bin_edges=np.histogram(a,bins=10,range=(0,1))
      adj_nb=sum(hist)
      hist=hist*1.0/sum(hist)*100.0

      spe_tag=spe[0]+"."+spe[10:13]
      # print spe_tag
      dict_spe_tag[spe]=spe_tag+"   ("+str(adj_nb)+")"
      dico_refine[spe]=hist

   fig, ax = plt.subplots()
   colors=['red','blue','yellow','purple','orange','green','brown','black','grey','coral','fuchsia','crimson','peru','olive','darkviolet','cyan']
   i=0
   for spe in dico_refine:
      lab=dict_spe_tag[spe]
      y=dico_refine[spe]
      x=np.arange(0.05,1.05,0.1)
      # print y
      # print x
      ax.plot(x,y,color=colors[i],label=lab)
      i=i+1

   handles, labels = ax.get_legend_handles_labels()
   hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
   handles2, labels2 = zip(*hl)

   # ax.legend(handles2,labels2)

   ax.set_xlabel("Score value interval")
   ax.set_ylabel("Frequency (%)")


   plt.xticks(np.arange(0, 1.1, .1))
   plt.yticks(np.arange(0, 81, 10))

   fig_name=OUTPUT+"/hist_BESST_scores_SPE_"+tag+".pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()








def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise



if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   file_scaff=open(argv[1],'r')
   output_graph_dir=argv[2]
   legend_tag=argv[3]


   mkdir_p(output_graph_dir)


   dict_species_list_vscore=dict()
   dict_species_list_dscore=dict()
   dict_species_list_mscore=dict()
   list_vscore=list()
   list_dscore=list()
   list_mscore=list()
   for line in file_scaff:

      r1=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', line)
      r2=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', line)
      species=""
      vscore=0.0
      dscore=0.0

      if r1:
         species=r1.group(1)
         if species!="#species":
            vscore=float(r1.group(7))
            dscore=float(r1.group(8))

            if not species in dict_species_list_vscore:
               dict_species_list_vscore[species]=list()
               dict_species_list_dscore[species]=list()
               dict_species_list_mscore[species]=list()

            list_vscore.append(vscore)
            list_dscore.append(dscore)
            list_mscore.append((vscore+dscore)/2.0)

            dict_species_list_vscore[species].append(vscore)
            dict_species_list_dscore[species].append(dscore)
            dict_species_list_mscore[species].append((vscore+dscore)/2.0)


      elif r2:
         species=r2.group(1)
         if species!="#species":
            vscore=float(r2.group(14))
            dscore=float(r2.group(15))

            if not species in dict_species_list_vscore:
               dict_species_list_vscore[species]=list()
               dict_species_list_dscore[species]=list()
               dict_species_list_mscore[species]=list()

            list_vscore.append(vscore)
            list_dscore.append(dscore)
            list_mscore.append((vscore+dscore)/2.0)

            dict_species_list_vscore[species].append(vscore)
            dict_species_list_dscore[species].append(dscore)
            dict_species_list_mscore[species].append((vscore+dscore)/2.0)

      else:
         exit("ERROR, the file is in a bad format!!!")


   hist_all(list_dscore,list_vscore,list_mscore,output_graph_dir,legend_tag)

   tag_vscore="vscore_"+legend_tag
   tag_dscore="dscore_"+legend_tag
   tag_mscore="mscore_"+legend_tag
   # hist_spe(dict_species_list_vscore,output_graph_dir,tag_vscore)
   # hist_spe(dict_species_list_dscore,output_graph_dir,tag_dscore)
   # hist_spe(dict_species_list_mscore,output_graph_dir,tag_mscore)

   hist_spe_norm(dict_species_list_vscore,output_graph_dir,tag_vscore)
   hist_spe_norm(dict_species_list_dscore,output_graph_dir,tag_dscore)
   hist_spe_norm(dict_species_list_mscore,output_graph_dir,tag_mscore)



   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))

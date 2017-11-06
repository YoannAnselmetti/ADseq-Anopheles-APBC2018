#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
###   Goal:
###      Plot distribution of gene families statistics
###
###   INPUT:
###      1- INPUT Gene Trees file
###         (data/INPUT_DATA/unrooted_raw_trees.nwk)
###         (data/GENE_TREES/unrooted_trees_filtered.nwk)
###         (datat/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk)
###      2- gene_TAG-species_name association file
###         (data/INPUT_DATA/name_geneID_18Anopheles)
###      3- OUTPUT_directory to store distribution graph of gene families
###         (figures/distrib_gene_trees/01-no_filter)
###         (figures/distrib_gene_trees/02-filter)
###         (figures/distrib_gene_trees/03-pNJ)
###      4- #Genes max to print in plot distribution
###         (40)
###
###   OUTPUT:
###      - Plots of gene families statistics
###
###   Name: distrib_gene_trees.py             Author: Yoann Anselmetti
###   Creation date: 2016/09/09               Last modification: 2017/01/27
###

from sys import argv
from re import search, sub, match
from os import close, rename, path, makedirs
from datetime import datetime
from shutil import rmtree
from Bio import Phylo   # Biopython library

import numpy as np 
from matplotlib import pyplot as plt

def file_len(fname):
   with open(fname) as f:
      for i, l in enumerate(f):
         pass
   return i+1


def replaceStringInFile(filePath,old_pattern,new_pattern):
   "replaces all string by a regex substitution"
   tempName = filePath+'~~~'
   inputFile = open(filePath)
   outputFile = open(tempName,'w')
   fContent = unicode(inputFile.read(), "utf-8")

   outText = sub(old_pattern, new_pattern, fContent)
   outputFile.write((outText.encode("utf-8")))

   outputFile.close()
   inputFile.close()

   rename(tempName, filePath)

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def cmpval(x,y):
    if x[1]>y[1]:
        return -1
    elif x[1]==y[1]:
        return 0
    else:
        return 1


# functions to plot distributions

def plot_species_on_gf_gene_distrib(SpeciesNb_geneNb,speciesNb_GFnb,OUTPUT):
   fig = plt.figure()
   ax1 = fig.add_subplot(111)

   list_GFnb=list()
   list_geneNb=list()
   i=1
   max_value=int(max(speciesNb_GFnb.keys()))
   while i<=max_value:
      if i in speciesNb_GFnb:
         list_GFnb.append(speciesNb_GFnb[i])
      else:
         list_GFnb.append(0)
      if i in SpeciesNb_geneNb:
         list_geneNb.append(SpeciesNb_geneNb[i])
      else:
         list_geneNb.append(0)
      i+=1

   ind = np.arange(max_value)
   width = 0.4

   ##################
   ### FOR axis 1 ###
   ##################
   p1 = ax1.bar(ind, list_GFnb, width, color='red')

   ax1.set_ylim(0,5000)
   ax1.yaxis.grid(True)
   ax1.set_yticks(np.arange(0,5001,500))
   ax1.set_ylabel("# Gene families",color="red")
   for tl in ax1.get_yticklabels():
      tl.set_color('red')

   # ax1.set_title("Distribution of #(Gene families) / #Species\nDistribution of #Genes in gene families with X species\n")
   ax1.set_xlim(-width,max_value+width)
   ax1.set_xticks(ind+width)
   ax1.set_xlabel("# Species")
   xTickMarks = [str(i) for i in range(1,max_value+1)]
   xtickNames = ax1.set_xticklabels(xTickMarks)

   ax2 = ax1.twinx()
   p1 = ax2.bar(ind+width, list_geneNb, width, color='blue')

   ax2.set_ylim(0,100000)
   ax2.yaxis.grid(True)
   ax2.set_yticks(np.arange(0,100001,10000))
   ax2.set_ylabel("# Genes",color="blue")
   for tl in ax2.get_yticklabels():
      tl.set_color('blue')


   plt.setp(xtickNames, rotation=67.5, fontsize=10)
   fig_name=OUTPUT+"/distrib_species_on_geneFamilies_or_gene.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def plot_species_gene_gf_distrib(geneNb_GFnb,speciesNb_GFnb,gene_max_plot,OUTPUT):
   fig = plt.figure()
   ax = fig.add_subplot(111)

   list_speNb=list()
   list_geneNb=list()
   i=1
   geneNb_in_GF_unprint=0
   GFnb_unprint=0
   max_value=int(max(geneNb_GFnb.keys()))
   while i<=max_value:
      if i in speciesNb_GFnb:
         list_speNb.append(speciesNb_GFnb[i])
      else:
         list_speNb.append(0)
      if i in geneNb_GFnb:
         list_geneNb.append(geneNb_GFnb[i])
         if i>gene_max_plot:
            geneNb_in_GF_unprint+=(geneNb_GFnb[i]*i)
            GFnb_unprint+=geneNb_GFnb[i]
      else:
         list_geneNb.append(0)
      i+=1

   ind = np.arange(max_value)
   width = 0.4
   ind2 = np.arange(gene_max_plot)

   p1 = ax.bar(ind, list_speNb, width, color='red')
   p2 = ax.bar(ind+width, list_geneNb, width, color='green')
   ax.legend((p1[0], p2[0]),('# Species', '# Genes'),loc='upper left', fontsize=15)

   ax.text(30.5, 3250, "Max gene family size: "+str(max_value), verticalalignment='center', horizontalalignment='center', color='green', fontsize=15)
   ax.text(30.5, 2250, "#Genes in the "+str(GFnb_unprint)+" gene families\nwith size >"+str(gene_max_plot)+": "+str(geneNb_in_GF_unprint), verticalalignment='center', horizontalalignment='center', color='green', fontsize=15)
   ax.set_xlim(-width,gene_max_plot+width)
   ax.set_ylim(0,5000)
   ax.yaxis.grid(True)
   ax.set_yticks(np.arange(0,5001,500))
   ax.set_xticks(ind2+width)
   ax.set_title("Distribution of #Gene families / #(Species or Genes)\n")
   ax.set_xlabel("# Species / # Genes")
   ax.set_ylabel("# Gene families")
   xTickMarks = [str(i) for i in range(1,gene_max_plot+1)]
   xtickNames = ax.set_xticklabels(xTickMarks)
   plt.setp(xtickNames, rotation=67.5, fontsize=10)
   fig_name=OUTPUT+"/distrib_species_or_gene_on_geneFamilies.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   INPUT_trees=argv[1]
   speciesName_geneID_file=open(argv[2],"r")
   OUTPUT_dir=argv[3]
   gene_max_plot=int(argv[4])

   sep="@"
   order_bool="prefix"

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_dir=path.normpath(OUTPUT_dir)

   # Remove OUTPUT_dir if existing
   if not path.exists(OUTPUT_dir):
      mkdir_p(OUTPUT_dir)

   log_file=open(OUTPUT_dir+"/log_file.txt","w")


   # Store association Gene-Species in "dict_gene_species"
   dict_geneID_speciesName={}
   for gene in speciesName_geneID_file:
      r=search("^([^\t ]*)[\t ]*([^\t ]*)\n$",gene)
      if r:
         speciesName=r.group(1)
         geneID=r.group(2)
         dict_geneID_speciesName[geneID]=speciesName
   speciesName_geneID_file.close()


#######################################################################################
### USE BIOPYTHON TO READ GENE TREES AND GET GENE PRESENT PRESENT IN EACH GENE TREE ###
#######################################################################################
   speciesNb_GFnb=dict()
   geneNb_GFnb=dict()
   SpeciesNb_geneNb=dict()
   spe_in_GFsize3=dict()
   speciesName_geneNb=dict()
   list_ASTEI=list()
   GFsize3=0
   i=0
   ASTEI_inc=0
   GF_onespe=0
   print "\nGene trees with size==1 (in species):"
   # Parse trees in gene trees file
   trees = Phylo.parse(INPUT_trees, 'newick')
   # For each tree
   for tree in trees:
      list_genes=list()
      list_species=list()
      i+=1
      # Give an ID to Gene families (Max: 9.999.999 gene families)
      GF_ID=""
      if i>=1000000 and i<=9999999:
         GF_ID="GF"+str(i)
      elif i>=100000:
         GF_ID="GF0"+str(i)
      elif i>=10000:
         GF_ID="GF00"+str(i)
      elif i>=1000:
         GF_ID="GF000"+str(i)
      elif i>=100:
         GF_ID="GF0000"+str(i)
      elif i>=10:
         GF_ID="GF00000"+str(i)
      elif i>=1:
         GF_ID="GF000000"+str(i)
      else:
         exit("ERROR too much gene families!!! (>=10.000.000)")


      # print GF_ID+":"
      # Get all clades of the gene tree to recover gene composing the gene tree
      for clade in tree.find_clades():
         if clade.name:
            # If clade correspond to a gene: add to the list of gene
            if match("^[A-Za-z]*[0-9]*$",clade.name) or match("^[A-Za-z]*[0-9]*-P[A-Z]$",clade.name):
               # print "\t"+clade.name
               for geneID in dict_geneID_speciesName:
                  if geneID in clade.name:
                     species=dict_geneID_speciesName[geneID]
                     if not species in speciesName_geneNb:
                        speciesName_geneNb[species]=0
                     speciesName_geneNb[species]+=1
                     list_genes.append(clade.name)
                     if not species in list_species:
                        list_species.append(species)
                     break
               # To analyze genes of Anopheles_stephensiI (ASTEI genes)
               if ("ASTEI" in clade.name):
                  if clade.name in list_ASTEI:
                     ASTEI_inc+=1
                     # print "Gene "+clade.name+" is present in several gene trees!!!"
                  else:
                     list_ASTEI.append(clade.name)

            elif match("^[A-Z][a-z]*_[a-z]*"+sep+"[A-Za-z]*[0-9]*$",clade.name) or match("^[A-Z][a-z]*_[a-z]*"+sep+"[A-Za-z]*[0-9]*-P[A-Z]$",clade.name) or match("^[A-Za-z]*[0-9]*"+sep+"[A-Z][a-z]*_[a-z]*$",clade.name) or match("^[A-Za-z]*[0-9]*-P[A-Z]"+sep+"[A-Z][a-z]*_[a-z]*$",clade.name):
               gene=""
               species=""
               if order_bool=="prefix":
                  species=str(clade.name).split(sep)[0]
                  gene=str(clade.name).split(sep)[1]
               elif order_bool=="postfix":
                  gene=str(clade.name).split(sep)[0]
                  species=str(clade.name).split(sep)[1]
               else:
                  exit("ERROR, parameter 5 should be equal to \"prefix\" or \"postfix\" !!!")

               # print "clade.name: "+clade.name
               # print "species: "+species
               # print "gene: "+gene

               if not species in speciesName_geneNb:
                  speciesName_geneNb[species]=0
               speciesName_geneNb[species]+=1
               list_genes.append(gene)
               if not species in list_species:
                  list_species.append(species)

               # To analyze genes of Anopheles_stephensiI (ASTEI genes)
               if ("ASTEI" in gene):
                  if gene in list_ASTEI:
                     ASTEI_inc+=1
                     # print "Gene "+clade.name+" is present in several gene trees!!!"
                  else:
                     list_ASTEI.append(gene)

      # To print gene of gene families of size 1 (in species
      if len(list_species)==1:
         GF_onespe+=1
         print "\t"+GF_ID+" => "+list_species[0]+":",
         for gene in list_genes:
            print "\t"+gene,
         print "\n",


      # Get number of species in current gene family and add it to speciesNb_GFnb dict()
      if not len(list_species) in speciesNb_GFnb:
         speciesNb_GFnb[len(list_species)]=1
      else:
         speciesNb_GFnb[len(list_species)]+=1

      # Get number of genes in current gene family and add it to geneNb_GFnb dict()
      if not len(list_genes) in geneNb_GFnb:
         geneNb_GFnb[len(list_genes)]=1
      else:
         geneNb_GFnb[len(list_genes)]+=1

      # Get number of genes in this gene family of size X species and add it to SpeciesNb_geneNb dict()
      if not len(list_species) in SpeciesNb_geneNb:
         SpeciesNb_geneNb[len(list_species)]=len(list_genes)
      else:
         SpeciesNb_geneNb[len(list_species)]+=len(list_genes)


      # Get occurence of species in gene families of size <=3 genes
      if len(list_genes)<=3:
         GFsize3+=1
         for spe in list_species:
            if not spe in spe_in_GFsize3:
               spe_in_GFsize3[spe]=1
            else:
               spe_in_GFsize3[spe]+=1


   print "\n=> There are "+str(i)+" gene families in file "+INPUT_trees
   print "\n=> There are "+str(GF_onespe)+" gene families with only 1 species!!!"
   print "\n=> There are "+str(ASTEI_inc)+" genes that are present in several trees in Anopheles_stephensiI!!!"
   log_file.write("=> There are "+str(i)+" gene families in file "+INPUT_trees+"\n")
   log_file.write("=> There are "+str(GF_onespe)+" gene families with only 1 species!!!\n")
   log_file.write("\n=> There are "+str(ASTEI_inc)+" genes that are present in several trees in Anopheles_stephensiI!!!\n")


   # Plot distribution graph
   plot_species_on_gf_gene_distrib(SpeciesNb_geneNb,speciesNb_GFnb,OUTPUT_dir)
   plot_species_gene_gf_distrib(geneNb_GFnb,speciesNb_GFnb,gene_max_plot,OUTPUT_dir)

   # Print table of distribution 
   i=1
   max_value=int(max(geneNb_GFnb.keys()))
   log_file.write("Value\t#Gene_families(/#species)\t#Gene_families(/#genes)\t#Genes(/gene_families_with_#species)\n")
   while i<=max_value:
      if i in speciesNb_GFnb:
         if i in geneNb_GFnb:
            log_file.write(str(i)+"\t"+str(speciesNb_GFnb[i])+"\t"+str(geneNb_GFnb[i])+"\t"+str(SpeciesNb_geneNb[i])+"\n")
         else:
            log_file.write(str(i)+"\t"+str(speciesNb_GFnb[i])+"\tNA\t"+str(SpeciesNb_geneNb[i])+"\n")
      else:
         if i in geneNb_GFnb:
            log_file.write(str(i)+"\tNA\t"+str(geneNb_GFnb[i])+"\tNA\n")
      i+=1


   print "\n=> Number of genes / species in gene families:"
   log_file.write("\nNumber of genes / species in gene families:\n")
   nb_gene_tot=0
   for spe in sorted(speciesName_geneNb):
      print "\t"+spe+"\t"+str(speciesName_geneNb[spe])
      log_file.write("\t"+spe+"\t"+str(speciesName_geneNb[spe])+"\n")
      nb_gene_tot+=speciesName_geneNb[spe]
   print "\tALL_species\t"+str(nb_gene_tot)
   log_file.write("\tALL_species\t"+str(nb_gene_tot)+"\n")


   print "\n=> Occurence of species present in gene families <= 3genes ("+str(GFsize3)+" gene families):"
   log_file.write("\nOccurence of species present in gene families <= 3genes ("+str(GFsize3)+" gene families):\n")
   L=spe_in_GFsize3.items()
   L.sort(cmpval)
   for spe in L:
      print "\t"+spe[0]+" => "+str(spe[1])
      log_file.write("\t"+spe[0]+" => "+str(spe[1])+"\n")


   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
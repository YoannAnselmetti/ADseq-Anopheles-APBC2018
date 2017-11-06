#!/usr/bin/env python
# -*- coding: utf-8 -*-
###                                                                          
###   Goal:                                                                  
###      To filter genes of species Anopheles_stephensiI taht contains genes with ID truncated and then present in several trees...
###
###   INPUT:                                                                 
###      1- INPUT Gene Trees file                                            
###         (data/INPUT_DATA/unrooted_trees_nofilt.nwk)
###      2- gene_TAG-species_name association file                           
###         (data/INPUT_DATA/name_geneID_18Anopheles)      
###      3- OUTPUT gene trees file                                           
###         (data/INPUT_DATA/unrooted_raw_trees.nwk) 
###      4- Minimum species number / gene tree                               
###         (1)                                                              
###                                                                          
###   OUTPUT:                                                                
###      - Filtered gene trees file with only gene form species present in parameter file 2                                 
###                                                                          
###   Name: filter_species_in_gene_trees.py       Author: Yoann Anselmetti
###   Creation date: 2016/09/09                   Last modification: 2017/10/26
###

from sys import argv
from re import search
from os import close, path, makedirs, remove
from datetime import datetime
from ete3 import Tree

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
   INPUT_trees=argv[1]
   speciesName_geneID_file=open(argv[2],"r")
   OUTPUT_trees_file=argv[3]
   species_threshold=int(argv[4])

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_dir=path.dirname(OUTPUT_trees_file)

   # Remove OUTPUT_dir if existing
   if not path.exists(OUTPUT_dir):
      mkdir_p(OUTPUT_dir)


   # Store association Gene-Species in "dict_gene_species"
   dict_geneID_speciesName={}
   for gene in speciesName_geneID_file:
      r=search("^([^\t ]*)[\t ]*([^\t ]*)\n$",gene)
      if r:
         speciesName=r.group(1)
         geneID=r.group(2)
         dict_geneID_speciesName[geneID]=speciesName
   speciesName_geneID_file.close()



  ######################
  ## with ETE3 module ##  => Function to prune gene trees!!!
  ######################

   input_file=open(INPUT_trees,"r")
   list_trees=list()
   # Browse gene trees file line by line
   for tree_line in input_file:
      # print tree_line
      bool_prune=False
      list_genes=list()
      list_species=list()
      tree_str=tree_line.replace("\n","")
      # print tree_str
      tree=Tree(tree_str)
      # Get list of extant genes in current gene tree
      for gene in tree.get_leaf_names():
         gene_present=False
         # Check if species of the gene is present in "speciesName_geneID_file". Else we have to prune current gene tree to filter this gene
         for geneID in dict_geneID_speciesName:
            if geneID in gene:
               gene_present=True
               list_genes.append(gene)
               if not dict_geneID_speciesName[geneID] in list_species:
                  list_species.append(dict_geneID_speciesName[geneID])
               break
         if not gene_present:
            bool_prune=True

      # Keep gene tree only if size (#genes) > 1
      if len(list_genes)>1:
         # Keep gene tree only if size (#species) > species_threshold
         if len(list_species)>=species_threshold:
            # If tree contains gene of a species that is not present in "speciesName_geneID_file" => Prune the tree
            if bool_prune:
               tree.prune(list_genes)
            # Write current gene tree in list_trees
            list_trees.append(tree.write(format=5))
   input_file.close()
   remove(INPUT_trees)

   output_file=open(OUTPUT_trees_file,"w")
   # Write list of gene trees in OUTPUT_trees_file
   for trees in list_trees:
      output_file.write(trees+"\n")


   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
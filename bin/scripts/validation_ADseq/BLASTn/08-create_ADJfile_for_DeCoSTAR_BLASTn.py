#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create Gene and Adjacency files fro ART-DeCo from Orthologous file with gene and exon position
###
###   INPUT:
###      1- INPUT file with gene position for all species
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/annotFile_Aalb_50pourc_allspecies)
###      2- INPUT scaffolding adjacencies file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/scaff_Aalb_50pourc_3_allspecies)
###      3- OUTPUT file path for Gene file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_albimanus/50pourc/INPUT_DeCoSTAR/adjacencies_Aara_50pourc)
###      4- Character separator between species name and gene ID
###         (@)
###      5- prefix/postfix boolean
###         (prefix or postfix)
###
###   OUTPUT:   (RUN in <1sec)
###      - INPUT files for ARt-DeCo_seq (Gene & Adjacency files)
###
###   Name: 08-create_ADJfile_for_DeCoSTAR_BLASTn.py  Author: Yoann Anselmetti
###   Creation date: 2016/04/01                       Last modification: 2016/11/17
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs


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
   SCAFF_file=argv[2]
   OUTPUT_Adj_ADseq=argv[3]
   OUTPUT_Adj_AD=argv[4]
   separator=argv[5]
   order_bool=argv[6]

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_Adj_ADseq))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   gene_file=open(GENE_file,'r')
   scaff_file=open(SCAFF_file,'r')
   output_adj_adseq=open(OUTPUT_Adj_ADseq,'w')
   output_adj_ad=open(OUTPUT_Adj_AD,'w')

   species=""
   contig=""
   gene_stored=""
   ori_stored="?"
   stop_stored=""
   gene2=""

   for gene in gene_file:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', gene)
      if r:
         name_species=r.group(1)
         contig_ID=r.group(2)
         gf_ID=r.group(3)
         gene_ID=r.group(4)
         orientation=r.group(5)
         start=r.group(6)
         stop=r.group(7)

         if name_species!="species":
            if name_species==species and contig_ID==contig:
               if order_bool=="prefix":
                  output_adj_adseq.write(name_species+separator+gene_stored+"\t"+name_species+separator+gene_ID+"\t"+ori_stored+"\t"+orientation+"\n")
                  output_adj_ad.write(name_species+separator+gene_stored+"\t"+name_species+separator+gene_ID+"\t"+ori_stored+"\t"+orientation+"\n")
               elif order_bool=="postfix":
                  output_adj_adseq.write(gene_stored+separator+name_species+"\t"+gene_ID+separator+name_species+"\t"+ori_stored+"\t"+orientation+"\n")
                  output_adj_ad.write(gene_stored+separator+name_species+"\t"+gene_ID+separator+name_species+"\t"+ori_stored+"\t"+orientation+"\n")
               else:
                  exit("ERROR, parameter 7 should be equal to \"postfix\" or \"prefix\" !!!")
               
               gene_stored=gene_ID
               ori_stored=orientation
               stop_stored=stop
            else:
               species=name_species
               contig=contig_ID
               gene_stored=gene_ID
               ori_stored=orientation
               stop_stored=stop
      else:
           exit("!!! ERROR line: "+gene+" is bad written !!!")
   gene_file.close()
   output_adj_ad.close()

   for scaff_adj in scaff_file:
      r=search('^([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\n]*)\n$', scaff_adj)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         ori_CTG1=r.group(4)
         ori_CTG2=r.group(5)
         gap=r.group(6)
         gf1=r.group(7)
         gf2=r.group(8)
         g1=r.group(9)
         g2=r.group(10)
         ori1=r.group(11)
         ori2=r.group(12)
         dist=r.group(13)
         link_score=r.group(14)
         disp_score=r.group(15)
         links_nb=r.group(16)

         if species!="species":
            score=(float(link_score)+float(disp_score))/2.0
            if order_bool=="prefix":
               output_adj_adseq.write(species+separator+g1+"\t"+species+separator+g2+"\t"+ori1+"\t"+ori2+"\t"+str(score)+"\n")
            elif order_bool=="postfix":
               output_adj_adseq.write(g1+separator+species+"\t"+g2+separator+species+"\t"+ori1+"\t"+ori2+"\t"+str(score)+"\n")
            else:
               exit("ERROR, parameter 7 should be equal to \"postfix\" or \"prefix\" !!!")
      else:
         exit("!!! ERROR line: "+scaff_adj+" is bad written !!!")
   scaff_file.close()
   output_adj_adseq.close()

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))

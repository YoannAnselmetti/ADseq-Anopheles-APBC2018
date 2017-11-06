#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Compare annotation gene files between original and after mapping of genes on genome assembly computed by minia
###
###   INPUT:
###      1- INPUT Original annotation file
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_albimanus/50pourc/annotG_Aalb_ORI_filt)
###      2- DeCo* file results after linearization 
###         (results/validation_ADseq/Aalb/50pourc/ADseq/validation_ADseq_Aalb_50pourc_Xtopo_Boltz_0.1.adjacencies.txt)
###      3- Species name
###         (Anopheles_albimanus)
###      4- Output file where stats will be stored for ADseq/ARt-DeCo
###         (results/validation_ADseq/stats/Aalb_50pourc_ADseq)
###      5- Orientation boolean to indicate if we have to take into account orientation in adjacencies or not
###         (T/t/True/true  or F/f/False/false)
###
###   OUTPUT:
###      - Give adjacencies that are in adjacencies file proposed by DeCo* (param 1) and are inconsistent with adajcencies in reference genome assembly (param 3)
###
###   Name: compare_newADJminia_ADJori_without_linearization.py      Author: Yoann Anselmetti
###   Creation date: 2016/07/18                                      Last modification: 2017/02/25
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
import errno

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


def rev_ori(ori_gene):
   if ori_gene=="-":
      return "+"
   elif ori_gene=="+":
      return "-"
   else:
      exit("ERROR, gene orientation is incorrect. Should be \"+\" or \"-\" and not \""+ori_gene+"\" !!!")


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   ORI_annotG_file=open(argv[1],'r')
   newAdjFILE=argv[2]
   select_spe=argv[3]
   output_stats=open(argv[4],'a')
   orientation=argv[5]

   bool_ori=False
   if orientation=="T" or orientation=="True" or orientation=="true" or orientation=="t" or orientation=="Y" or orientation=="Yes" or orientation=="yes" or orientation=="y":
      bool_ori=True
   elif orientation=="F" or orientation=="False" or orientation=="false" or orientation=="f" or orientation=="N" or orientation=="No" or orientation=="no" or orientation=="n":
      bool_ori=False
   else:
      exit("ERROR, parameter 4 for orientation boolean should be equal to \"true\" or \"false\" and NOT: "+orientation+" !!!")



   ADJ=namedtuple("ADJ",["g1","g2","ori1","ori2"])
   threshold_values=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99]

######################
### Determine list of adjacencies in reference genome assembly of selected species
######################
   gene_adj=dict()
   species=""
   contig=""
   gene_store=""
   adj_tot_ori=0
   for gene in ORI_annotG_file:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', gene)
      if r:
         name_species=r.group(1)
         contig_ID=r.group(2)
         gf_ID=r.group(3)
         gene_ID=r.group(4)
         gene_ori=r.group(5)
         gene_start=r.group(6)
         gene_end=r.group(7)


         spe_gene=name_species+"@"+gene_ID

         if name_species==select_spe:
            gene_adj[spe_gene]=list()
            if name_species==species and contig_ID==contig:
               adj_tot_ori+=1

               adj=ADJ(gene_store,spe_gene,ori_store,gene_ori)
               gene_adj[gene_store].append(adj)
               gene_adj[spe_gene].append(adj)

               gene_store=spe_gene
               ori_store=gene_ori
            else:
               species=name_species
               contig=contig_ID
               gene_store=spe_gene
               ori_store=gene_ori
      else:
           exit("ERROR line: "+gene+" is bad written")
   ORI_annotG_file.close()



   for threshold in threshold_values:
   ######################
   ### Analysed new extant adjacencies predicted by DeCo* in MINIA genome assembly of selected species
   ######################
      TP=0
      TP_ori=0
      FP=0
      CFP=0
      new_adj=0
      obs_adj=0
      incorrect_ori=0
      newAdj_file=open(newAdjFILE,'r')
      for adj in newAdj_file:
         r=search('^([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ \n]*)\n$', adj)
         if r:
            species=r.group(1)
            gene1=r.group(2)
            gene2=r.group(3)
            ori1=r.group(4)
            ori2=r.group(5)
            prior_score=r.group(6)
            post_score=r.group(7)

            species_name=gene1.split("@")[0]

            # cur_adj=ADJ(gene1,gene2,ori1,ori2)
            if species_name==select_spe:
               # print adj
               # If prior score==1, then ADJ is an observed one
               if float(prior_score)==1.0:
                  obs_adj+=1
                  # Analyze if the adajcency is well present in reference assembly
                  present=False
                  correc_ori=False
                  for adj in gene_adj[gene1]:
                     g1=adj.g1
                     g2=adj.g2
                     o1=adj.ori1
                     o2=adj.ori2
                     if gene1==g1:
                        if g2==gene2:
                           present=True
                           if ori1==o1 and ori2==o2:
                              correc_ori=True
                           break
                     elif gene1==g2:
                        if g1==gene2:
                           present=True
                           o1=rev_ori(o1)
                           o2=rev_ori(o2)
                           if ori1==o2 and ori2==o1:
                              correc_ori=True
                           break
                  if not present:
                     exit("ERROR: ADJ "+gene1+"-"+gene2+" IS NOT PRESENT IN REFERENCE GENOME ASSEMBLY !!!")
                  else:
                     if not correc_ori:
                        incorrect_ori+=1
                        print "Genes "+gene1+" and "+gene2+" are incorrectly orientated !!!"
               # If score prior!=1, then ADJ is new adajcencies (scaff or not)
               else:
                  if float(post_score)>0.0:
                     if float(post_score)>=threshold:
                        new_adj+=1
                        # Analyze if the adajcency is well present in reference assembly
                        present=False
                        correc_ori=False
                        for adj in gene_adj[gene1]:
                           g1=adj.g1
                           g2=adj.g2
                           o1=adj.ori1
                           o2=adj.ori2
                           if gene1==g1:
                              if g2==gene2:
                                 present=True
                                 if ori1==o1 and ori2==o2:
                                    correc_ori=True
                                 break

                           elif gene1==g2:
                              if g1==gene2:
                                 present=True
                                 o1=rev_ori(o1)
                                 o2=rev_ori(o2)
                                 if ori1==o2 and ori2==o1:
                                    correc_ori=True
                                 break



                        if present:
                           if bool_ori:
                              if correc_ori:
                                 TP+=1
                              else:
                                 FP+=1
                                 CFP+=1
                           else:
                              TP+=1
                        else:
                           FP+=1
                           # print "NEW ADJ: "+gene1+" "+gene2+" -> NOT PRESENT",
                           if (len(gene_adj[gene1])==2 or len(gene_adj[gene2])==2):
                              CFP+=1
                              #    print "AND INCONSISTENT WITH ADJ",
                           # print "IN INITIAL GENOME ASSEMBLY"

      newAdj_file.close()


      adj_to_find=adj_tot_ori-obs_adj
      print "\nSTATS on validation of DeCo* for genome fragmentation simulation of species "+select_spe+":"
      print "\t- "+str(adj_to_find)+" observed adjacencies in file reference genome assembly NOT present in MINIA genome assembly before DeCo* execution for species "+select_spe+"\n"
      if new_adj==0:
         print "THERE IS NO NEW ADJACENCIES FOR THIS EXPERIMENT"
         output_stats.write(threshold+"\t0\t"+str(adj_to_find)+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
      else:
         print "\t- On the "+str(new_adj)+" new adjacencies for species "+select_spe+":"
         print "\t\t+ "+str(TP)+" were present in file: "+argv[1] 
         print "\t\t+ "+str(FP)+" were NOT present in file: "+argv[1]
         if FP!=0:
            print "\t\t  => "+str(CFP)+" are inconsistent with ADJ present in "+argv[1]

            CFP_FP=0.0
            FN=adj_to_find-TP

            FFP=FP-CFP

            precision="%.2f" % float(float(TP)/float(new_adj)*100.0)
            recall="%.2f" % float(float(TP)/float(adj_to_find)*100.0)

            precision_CFP="%.2f" % float(float(TP+FFP)/float(new_adj)*100.0)
            recall_CFP="%.2f" % float(float(TP+FFP)/float(adj_to_find)*100.0)


            print "\nPrecision= "+str(precision)+"%"
            if FP!=0:
               CFP_FP="%.2f" % float(float(CFP)/float(FP)*100.0)
               print "ratio CFP/FP= "+str(CFP_FP)+"%"
            print "Recall= "+str(recall)+"%"


            if precision==100.0:
               output_stats.write(str(threshold)+"\t"+str(new_adj)+"\t"+str(adj_to_find)+"\t"+str(FN)+"\t"+str(TP)+"\t"+str(FP)+"\t"+str(CFP)+"\tNA\t"+str(recall)+"\t"+str(precision)+"\t"+str(recall_CFP)+"\t"+str(precision_CFP)+"\n")
            else:
               output_stats.write(str(threshold)+"\t"+str(new_adj)+"\t"+str(adj_to_find)+"\t"+str(FN)+"\t"+str(TP)+"\t"+str(FP)+"\t"+str(CFP)+"\t"+str(CFP_FP)+"\t"+str(recall)+"\t"+str(precision)+"\t"+str(recall_CFP)+"\t"+str(precision_CFP)+"\n")

      end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
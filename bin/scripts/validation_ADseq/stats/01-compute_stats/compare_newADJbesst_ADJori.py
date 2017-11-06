#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Compare annotation gene files between original and after mapping of genes on genome assembly computed by minia
###
###   INPUT:
###      1- INPUT Original annotation file
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_albimanus/50pourc/annotG_Aalb_ORI_filt)
###      2- DeCo* file results without linearization 
###         (results/validation_ADseq/Aalb/50pourc/ADseq/validation_ADseq_Aalb_50pourc_Xtopo_Boltz_0.1.adjacencies.txt)
###      3- Species name
###         (Anopheles_albimanus)
###      4- AGP file BESST results
###         (data/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k50/50pourc/Anopeheles_albimanus/BESST_output/pass1/info-pass1.agp)
###      5- BESST new adjacencies results file 
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_albimanus/50pourc/scaff_Aalb_50pourc_3)
###      6- Species tag
###         (Aalb)
###      7- Sequencing tag
###         (50pourc)
###      8- OUTPUT directory where all stats will be stored
###         (results/validation_ADseq/stats/)
###      9- Orientation boolean to indicate if we have to take into account orientation in adjacencies or not
###         (T/t/True/true  or F/f/False/false)
###
###   OUTPUT:
###      - Give adjacencies that are in adjacencies file proposed by DeCo* (param 1) and are inconsistent with adajcencies in reference genome assembly (param 3)
###
###   Name: compare_newADJbesst_ADJori.py      Author: Yoann Anselmetti
###   Creation date: 2017/01/05                Last modification: 2017/02/25
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


def write_listADJ(listADJ,outFile,bool_ori):
   out=open(outFile,"w")
   for adj in listADJ:
      g1=adj.g1
      g2=adj.g2
      ori1=adj.ori1
      ori2=adj.ori2
      if bool_ori:
         out.write(g1+"\t"+g2+"\t"+ori1+"\t"+ori2+"\n")
      else:
         out.write(g1+"\t"+g2+"\n")

################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   ORI_annotG_file=open(argv[1],'r')
   newAdj_file=open(argv[2],'r')
   select_spe=argv[3]
   agp_file=open(argv[4],'r')
   BESST_FILE=argv[5]
   spe_tag=argv[6]
   seq_tag=argv[7]
   OUTPUT=argv[8]
   orientation=argv[9]

   bool_ori=False
   if orientation=="T" or orientation=="True" or orientation=="true" or orientation=="t" or orientation=="Y" or orientation=="Yes" or orientation=="yes" or orientation=="y":
      bool_ori=True
   elif orientation=="F" or orientation=="False" or orientation=="false" or orientation=="f" or orientation=="N" or orientation=="No" or orientation=="no" or orientation=="n":
      bool_ori=False
   else:
      exit("ERROR, parameter 4 for orientation boolean should be equal to \"true\" or \"false\" and NOT: "+orientation+" !!!")


   mkdir_p(OUTPUT)

   ADJ=namedtuple("ADJ",["g1","g2","ori1","ori2"])
   threshold_values=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99]

######################
### Determine list of adjacencies in reference genome assembly of selected species
######################
   list_FN=list()
   gene_adj=dict()
   species=""
   contig=""
   gene_store=""
   gene2=""
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

               list_FN.append(adj)

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



######################
### Get observed adjacencies in minia assembly to determine the number of adjacencies to find (fragmented adjacencies)
######################
   obs_adj=0
   obs_adj_supp=0
   incorrect_ori=0
   for adj in newAdj_file:
      r=search('^([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*)[ \t]([^ \t\n]*)\n$', adj)
      if r:
         species=r.group(1)
         gene1=r.group(2)
         gene2=r.group(3)
         ori1=r.group(4)
         ori2=r.group(5)
         prior_score=r.group(6)
         post_score=r.group(7)
         scaff=r.group(8)

         species_name=gene1.split("@")[0]

         # cur_adj=ADJ(gene1,gene2,ori1,ori2)
         if species_name==select_spe:
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

                        list_FN.remove(adj)
                        obs_adj_supp+=1

                        break

                  elif gene1==g2:
                     if g1==gene2:
                        present=True
                        o1=rev_ori(o1)
                        o2=rev_ori(o2)
                        if ori1==o2 and ori2==o1:
                           correc_ori=True

                        list_FN.remove(adj)
                        obs_adj_supp+=1

                        break

               if not present:
                  exit("ERROR: ADJ "+spe_gene1+"-"+spe_gene2+" IS NOT PRESENT IN REFERENCE GENOME ASSEMBLY !!!")
               else:
                  if not correc_ori:
                     incorrect_ori+=1
                     print "Genes "+gene1+" and "+gene2+" are incorrectly orientated !!!"
   newAdj_file.close()



   list_TP=list()
   list_FP_CFP=list()
   list_CFP=list()
   list_newADJ=list()
   for adj_line in agp_file:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', adj_line)
      if r:
         species=r.group(1)
         gene1=r.group(2)
         gene2=r.group(3)
         ori1=r.group(4)
         ori2=r.group(5)

         if species==select_spe:
            present=False
            correc_ori=False
            spe_gene1=species+"@"+gene1
            spe_gene2=species+"@"+gene2

            cur_adj=ADJ(spe_gene1,spe_gene2,ori1,ori2)
            list_newADJ.append(cur_adj)

            for adj in gene_adj[spe_gene1]:
               g1=adj.g1
               g2=adj.g2
               o1=adj.ori1
               o2=adj.ori2
               if spe_gene1==g1:
                  if g2==spe_gene2:
                     present=True
                     if ori1==o1 and ori2==o2:
                        correc_ori=True
                     break

               elif spe_gene1==g2:
                  if g1==spe_gene2:
                     present=True
                     o1=rev_ori(o1)
                     o2=rev_ori(o2)
                     if ori1==o2 and ori2==o1:
                        correc_ori=True
                     break

            if present:
               if bool_ori:
                  if correc_ori:
                     if adj in list_FN:
                        list_FN.remove(adj)
                        list_TP.append(cur_adj)
                     # else:
                     #    print cur_adj
                     #    print "\t=>",
                     #    print adj
                  else:
                     list_CFP.append(cur_adj)


               else:
                  if adj in list_FN:
                     list_TP.append(cur_adj)
                     list_FN.remove(adj)
                  # else:
                  #    print cur_adj
                  #    print "\t=>",
                  #    print adj
                  
            else:
               # print "NEW ADJ: "+spe_gene1+" "+spe_gene2+" -> NOT PRESENT",
               if (len(gene_adj[spe_gene1])==2 or len(gene_adj[spe_gene2])==2):
                  list_CFP.append(cur_adj)
                  #    print "AND INCONSISTENT WITH ADJ",
               else:
                  list_FP_CFP.append(cur_adj)
               # print "IN INITIAL GENOME ASSEMBLY"
               
   agp_file.close()


   mkdir_p(OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/TP")
   mkdir_p(OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/FP-CFP")
   mkdir_p(OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/CFP")
   mkdir_p(OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/ALL")
   mkdir_p(OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/FN")

   TP_file=OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/TP/BESST_complete"
   FP_CFP_file=OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/FP-CFP/BESST_complete"
   CFP_file=OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/CFP/BESST_complete"
   ALL_file=OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/ALL/BESST_complete"
   FN_file=OUTPUT+"/Venn_diagram/"+spe_tag+"/"+seq_tag+"/FN/BESST_complete"

   write_listADJ(list_TP,TP_file,bool_ori)
   write_listADJ(list_FP_CFP,FP_CFP_file,bool_ori)
   write_listADJ(list_CFP,CFP_file,bool_ori)
   write_listADJ(list_newADJ,ALL_file,bool_ori)
   write_listADJ(list_FN,FN_file,bool_ori)


   TP_besst=len(list_TP)
   FP_besst=len(list_CFP)+len(list_FP_CFP)
   CFP_besst=len(list_CFP)
   new_adj_besst=len(list_newADJ)


   new_adj_besst=len(list_newADJ)
####################################
### PRINT STATS PRECISION/RECALL ###
####################################
   output_stats=open(OUTPUT+"/"+spe_tag+"_"+seq_tag+"_BESST",'a')
   adj_to_find=adj_tot_ori-obs_adj
   print "\n\nFOR BESST (complete run)"
   if new_adj_besst==0:
      print "THERE IS NO NEW ADJACENCIES FOR THIS EXPERIMENT"
      output_stats.write("COMPLET\t0\t"+str(adj_to_find)+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
   else:
      print "\t- On the "+str(new_adj_besst)+" new adjacencies for species "+select_spe+":"
      print "\t\t+ "+str(TP_besst)+" were present in file: "+argv[1]
      print "\t\t+ "+str(FP_besst)+" were NOT present in file: "+argv[1]
      if FP_besst!=0:
         print "\t\t  => "+str(CFP_besst)+" are inconsistent with ADJ present in "+argv[1]

         CFP_FP_besst=0.0
         FN_besst=adj_to_find-TP_besst

         print "FN_besst= "+str(FN_besst)+" == len(list_FN)= "+str(len(list_FN))

         FFP_besst=FP_besst-CFP_besst

         precision_besst="%.2f" % float(float(TP_besst)/float(new_adj_besst)*100.0)
         recall_besst="%.2f" % float(float(TP_besst)/float(adj_to_find)*100.0)


         precision_besst_CFP="%.2f" % float(float(TP_besst+FFP_besst)/float(new_adj_besst)*100.0)
         recall_besst_CFP="%.2f" % float(float(TP_besst+FFP_besst)/float(adj_to_find)*100.0)


         print "\nPrecision= "+str(precision_besst)+"%"
         if FP_besst!=0:
            CFP_FP_besst="%.2f" % float(float(CFP_besst)/float(FP_besst)*100.0)
            print "ratio CFP/FP= "+str(CFP_FP_besst)+"%"
         print "Recall= "+str(recall_besst)+"%"

         if precision_besst==100.0:
            output_stats.write("COMPLET\t"+str(new_adj_besst)+"\t"+str(adj_to_find)+"\t"+str(FN_besst)+"\t"+str(TP_besst)+"\t"+str(FP_besst)+"\t"+str(CFP_besst)+"\tNA\t"+str(recall_besst)+"\t"+str(precision_besst)+"\t"+str(recall_besst_CFP)+"\t"+str(precision_besst_CFP)+"\n")
         else:
            output_stats.write("COMPLET\t"+str(new_adj_besst)+"\t"+str(adj_to_find)+"\t"+str(FN_besst)+"\t"+str(TP_besst)+"\t"+str(FP_besst)+"\t"+str(CFP_besst)+"\t"+str(CFP_FP_besst)+"\t"+str(recall_besst)+"\t"+str(precision_besst)+"\t"+str(recall_besst_CFP)+"\t"+str(precision_besst_CFP)+"\n")


   for threshold in threshold_values:

######################
### Analysed new extant adjacencies proposed by BESST in MINIA genome assembly of selected species
######################
      TP_besst=0
      FP_besst=0
      CFP_besst=0
      new_adj_besst=0
      BESST_file=open(BESST_FILE,'r')
      for adj_line in BESST_file:
         r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)[\t]([^\t\n]*)\n$', adj_line)
         if r:
            species=r.group(1)
            ctg1=r.group(2)
            ctg2=r.group(3)
            oriCTG1=r.group(4)
            oriCTG2=r.group(5)
            distCTG=r.group(6)
            gf1=r.group(7)
            gf2=r.group(8)
            gene1=r.group(9)
            gene2=r.group(10)
            ori1=r.group(11)
            ori2=r.group(12)
            distG=r.group(13)
            vscore=r.group(14)
            dscore=r.group(15)
            links=r.group(16)


            if species==select_spe:
               if (float(vscore)+float(dscore)/2.0)>=threshold:
                  new_adj_besst+=1
                  present=False
                  correc_ori=False
                  spe_gene1=species+"@"+gene1
                  spe_gene2=species+"@"+gene2

                  for adj in gene_adj[spe_gene1]:
                     g1=adj.g1
                     g2=adj.g2
                     o1=adj.ori1
                     o2=adj.ori2
                     if spe_gene1==g1:
                        if g2==spe_gene2:
                           present=True
                           if ori1==o1 and ori2==o2:
                              correc_ori=True
                           break

                     elif spe_gene1==g2:
                        if g1==spe_gene2:
                           present=True
                           o1=rev_ori(o1)
                           o2=rev_ori(o2)
                           if ori1==o2 and ori2==o1:
                              correc_ori=True
                           break

                  if present:
                     if bool_ori:
                        if correc_ori:
                           TP_besst+=1
                        else:
                           FP_besst+=1
                           CFP_besst+=1
                     else:
                        TP_besst+=1
                  else:
                     FP_besst+=1
                     # print "NEW ADJ: "+spe_gene1+" "+spe_gene2+" -> NOT PRESENT",
                     if (len(gene_adj[spe_gene1])==2 or len(gene_adj[spe_gene2])==2):
                        CFP_besst+=1
                        #    print "AND INCONSISTENT WITH ADJ",
                     # print "IN INITIAL GENOME ASSEMBLY"

      BESST_file.close()


####################################
### PRINT STATS PRECISION/RECALL ###
####################################
      adj_to_find=adj_tot_ori-obs_adj
      print "\n\nFOR BESST (with support >"+str(threshold)+"):"
      if new_adj_besst==0:
         print "THERE IS NO NEW ADJACENCIES FOR THIS EXPERIMENT"
         output_stats.write(str(threshold)+"\t0\t"+str(adj_to_find)+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
      else:
         print "\t- On the "+str(new_adj_besst)+" new adjacencies for species "+select_spe+":"
         print "\t\t+ "+str(TP_besst)+" were present in file: "+argv[1]
         print "\t\t+ "+str(FP_besst)+" were NOT present in file: "+argv[1]
         if FP_besst!=0:
            print "\t\t  => "+str(CFP_besst)+" are inconsistent with ADJ present in "+argv[1]

         CFP_FP_besst=0.0
         FN_besst=adj_to_find-TP_besst

         FFP_besst=FP_besst-CFP_besst

         precision_besst="%.2f" % float(float(TP_besst)/float(new_adj_besst)*100.0)
         recall_besst="%.2f" % float(float(TP_besst)/float(adj_to_find)*100.0)


         precision_besst_CFP="%.2f" % float(float(TP_besst+FFP_besst)/float(new_adj_besst)*100.0)
         recall_besst_CFP="%.2f" % float(float(TP_besst+FFP_besst)/float(adj_to_find)*100.0)


         print "\nPrecision= "+str(precision_besst)+"%"
         if FP_besst!=0:
            CFP_FP_besst="%.2f" % float(float(CFP_besst)/float(FP_besst)*100.0)
            print "ratio CFP/FP= "+str(CFP_FP_besst)+"%"
         print "Recall= "+str(recall_besst)+"%"

         if precision_besst==100.0:
            output_stats.write(str(threshold)+"\t"+str(new_adj_besst)+"\t"+str(adj_to_find)+"\t"+str(FN_besst)+"\t"+str(TP_besst)+"\t"+str(FP_besst)+"\t"+str(CFP_besst)+"\tNA\t"+str(recall_besst)+"\t"+str(precision_besst)+"\t"+str(recall_besst_CFP)+"\t"+str(precision_besst_CFP)+"\n")
         else:
            output_stats.write(str(threshold)+"\t"+str(new_adj_besst)+"\t"+str(adj_to_find)+"\t"+str(FN_besst)+"\t"+str(TP_besst)+"\t"+str(FP_besst)+"\t"+str(CFP_besst)+"\t"+str(CFP_FP_besst)+"\t"+str(recall_besst)+"\t"+str(precision_besst)+"\t"+str(recall_besst_CFP)+"\t"+str(precision_besst_CFP)+"\n")


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
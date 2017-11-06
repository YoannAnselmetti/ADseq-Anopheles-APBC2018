#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create scaffolding graph file for creation of instances for AD/ADseq
###   INPUT:
###      1- Pre-scaffolding adjacencies file
###         (Validation_ADseq/coverage/ADseq/preScaff_ORI_3)
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/preScaff_Aara_50pourc_3_modif)
###      2- CTG file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/CTG_Aara_ORI)
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/CTG_Aara_50pourc_filt)
###      3- add gene to CTG scaffolds file
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/scaff_Aara_ORI_3)
###         (Validation_ADseq/coverage/ADseq/BLASTn/Anopheles_arabiensis/50pourc/scaff_Aara_50pourc_3)
###
###   OUTPUT:   (RUN in ~1h for ORI & ~10s for others)
###      - Create scaffolding adj file for creation of Instances for ADseq
###        => Take only 1 solution for each couple of gene adjacencies
###           (Same if different orientations possible)
###   Name: 06c-create_scaff_adj_file_final_BLASTn.py  Author: Yoann Anselmetti
###   Creation date: 2015/11/06                        Last modification: 2017/03/15
###

from sys import argv, stdout
from re import search
from os import close
from collections import namedtuple   #New in version 2.6
from datetime import datetime





def get_gene_infos(ori,dist,ctg,dict_ID_ctg,contigEXT_file,preSCAFF_file):
   GF_gene=""
   gene=""
   ori_gene=""
   # if contig ctg2 in "+" orientation: gene=g1
   if (ori=="+"):
      GF_gene=dict_ID_ctg[ctg].gf1
      gene=dict_ID_ctg[ctg].g1
      ori_gene=dict_ID_ctg[ctg].ori1
      dist+=dict_ID_ctg[ctg].start
   # if contg ctg in "-" orientation: gene=g2
   elif (ori=="-"):
      GF_gene=dict_ID_ctg[ctg].gf2
      gene=dict_ID_ctg[ctg].g2
      dist+=dict_ID_ctg[ctg].size-dict_ID_ctg[ctg].end
      # Get correct orientation of gene2
      oriG=dict_ID_ctg[ctg].ori2
      if (oriG=="-"):
         ori_gene="+"
      elif (oriG=="+"):
         ori_gene="-"
      elif (oriG=="?"):
         ori_gene="?"
      else:
         exit("ERROR on gene orientation in contigs extremities file "+contigEXT_file+" (Should be \"+\" or \"-\")")
   else:
      exit("ERROR on gene orientation in scaffolding graph file "+preSCAFF_file+" (Should be \"+\" or \"-\")")

   return GF_gene,gene,ori_gene,dist



if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   preSCAFF_file=argv[1]
   contigEXT_file=argv[2]
   OUTPUT_scaff_file=argv[3]

   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   ADJ=namedtuple("ADJ",["spe","g1","g2"])
   ADJORI=namedtuple("ADJORI",["spe","g1","g2","ori1","ori2"])
   CTG=namedtuple("CTG",["spe","ctg","size","gf1","g1","ori1","start","gf2","g2","ori2","end"])
   EDGE=namedtuple("EDGE",["spe","ctg1","ctg2","oriC1","oriC2","gap","gf1","gf2","g1","g2","oriG1","oriG2","dist","vscore","dscore","link"])

#######
### STORE CTG INFOS IN dict_spe_ID_ctg
#######
   print "1/ Store CTG infos.. ",
   dict_spe_ID_ctg={}
   # When get a contigs pairs (edge scaffolding link) => Get genes that are linked by scaffolding graph with the distance
   with open(contigEXT_file,'r') as contig_file:
      for line in contig_file:
         spe=""
         contig=""
         contig_size=""
         GF1=""
         g1=""
         oriG1=""
         start_g1=""
         GF2=""
         g2=""
         oriG2=""
         stop_g2=""
         r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         r2=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r:
            spe=r.group(1)
            contig=r.group(2)
            contig_size=r.group(3)
            GF1=r.group(4)
            g1=r.group(5)
            oriG1=r.group(6)
            start_g1=r.group(7)
            GF2=r.group(8)
            g2=r.group(9)
            oriG2=r.group(10)
            stop_g2=r.group(11)
         elif r2:
            spe=r2.group(1)
            contig=r2.group(2)
            contig_size=r2.group(3)
            contig_geneNb=r2.group(4)
            GF1=r2.group(5)
            g1=r2.group(6)
            oriG1=r2.group(7)
            start_g1=r2.group(8)
            GF2=r2.group(9)
            g2=r2.group(10)
            oriG2=r2.group(11)
            stop_g2=r2.group(12)
         else:
            exit("ERROR in line "+line+" of file "+contigEXT_file+" !!!")

         if spe!="#species":
            if contig_size=="?":
               print "\n\t=> Contig "+contig+" is not present in FASTA file assembly of species "+spe
            else:
               ctg=CTG(spe,contig,int(contig_size),GF1,g1,oriG1,int(start_g1),GF2,g2,oriG2,int(stop_g2))
               if not spe in dict_spe_ID_ctg:
                  dict_spe_ID_ctg[spe]=dict()
               dict_spe_ID_ctg[spe][contig]=ctg

   contig_file.close()
   print "DONE"

   # for spe in dict_spe_ID_ctg:
   #    print spe+":"
   #    for ctg in dict_spe_ID_ctg[spe]:
   #       print "\t"+ctg+":",
   #       print dict_spe_ID_ctg[spe][ctg]


#######
### CREATE SCAFFOLDING GRAPH FILE WITH GENES INVOLVED IN SCAFFFOLDING EDGES INSTEAD CONTIGS
#######
   print "2/ Store all scaffolding gene adjacencies present in file "+preSCAFF_file+"...",
   stdout.flush()
   dict_spe_edge_scaff=dict()
   Nb_scaff_edge_tot=0
   Nb_scaff_edge_kept=0
   # Browse edge from scaffolding graph file and store them in "dict_spe_edge_scaff[spe]" (1 scaffolding graph file for ALL species)
   with open(preSCAFF_file,'r') as pre_scaff_file:
      for line in pre_scaff_file:
         r1=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r1:
            species=r1.group(1)
            ctg1=r1.group(2)
            ctg2=r1.group(3)
            ori1=r1.group(4)
            ori2=r1.group(5)
            gap=r1.group(6)
            vscore=r1.group(7)
            dscore=r1.group(8)
            links_nb=r1.group(9)

            if species!="species":
               Nb_scaff_edge_tot+=1
               # distance between gene1 and gene2 => If gap is a negative distance == Size of the overlap? -> distance = gap + d(g1-ov) + d(ov-g2)
               dist=float(gap)
               bool_ctg1=False
               bool_ctg2=False

               GF_gene1=""
               gene1=""
               ori_gene1=""
               GF_gene2=""
               gene2=""
               ori_gene2=""
               # If ctg1 is present in CTG_file: Get information on gene involved in the linked between ctg1 and ctg2
               if ctg1 in dict_spe_ID_ctg[species]:
                  bool_ctg1=True
                  # Get infos for gene1
                  GF_gene1,gene1,ori_gene1,dist = get_gene_infos(ori1,dist,ctg1,dict_spe_ID_ctg[species],contigEXT_file,preSCAFF_file)

               # If ctg2 is present in contigEXT_file: Get information on gene involved in the linked between ctg1 and ctg2
               if ctg2 in dict_spe_ID_ctg[species]:
                  bool_ctg2=True
                  # Get infos for gene2
                  GF_gene2,gene2,ori_gene2,dist = get_gene_infos(ori2,dist,ctg2,dict_spe_ID_ctg[species],contigEXT_file,preSCAFF_file)

               score1='{0:.12f}'.format(float(vscore))
               score2='{0:.12f}'.format(float(dscore))


               # If the 2 contigs are present in file contigEXT_file: Print scaffolding adj in OUTPUT_scaff_file
               if bool_ctg1 and bool_ctg2:
                  edge=EDGE(species,ctg1,ctg2,ori1,ori2,float(gap),GF_gene1,GF_gene2,gene1,gene2,ori_gene1,ori_gene2,float(dist),score1,score2,int(links_nb))
                  if not species in dict_spe_edge_scaff:
                     dict_spe_edge_scaff[species]=dict()
                  adj=ADJORI(species,gene1,gene2,ori_gene1,ori_gene2)

                  # NEED TO ORDER GENES IN ADJ IF WE DON'T CONSIDER GENE ORIENTATION IN ADJ!!!
                  # adj=ADJ(species,gene1,gene2)
                  if adj in dict_spe_edge_scaff[species]:
                     print ""
                     print adj,
                     print "=> present in several copies (CTG1: "+ctg1+" & CTG2: "+ctg2+")"
                     exit()

                  dict_spe_edge_scaff[species][adj]=edge
 
         else:
            exit("ERROR in line:\n\t"+line+"of file "+preSCAFF_file+" !!!")
   pre_scaff_file.close()
   print "DONE"


#########################################################
### WRITE OUTPUT_FILE => FINAL SCAFFOLDING GRAPH FILE ###
#########################################################
   print "3/ Write final scaffolding gene adjacencies file... ",
   stdout.flush()
   nb_scaff_kept=0
   scaff_gene_file=open(OUTPUT_scaff_file,"w")
   scaff_gene_file.write("species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tctg1-ctg2_dist\tgene1_family\tgene2_family\tgene1\tgene2\torientation_gene1\torientation_gene2\tgene1-gene2_dist\tvscore\tdscore\t#links\n")
   for spe in sorted(dict_spe_edge_scaff):
      nb_scaff_kept+=len(dict_spe_edge_scaff[spe])
      for adj in sorted(dict_spe_edge_scaff[spe]):
         scaff_gene_file.write(dict_spe_edge_scaff[spe][adj].spe+"\t"+dict_spe_edge_scaff[spe][adj].ctg1+"\t"+dict_spe_edge_scaff[spe][adj].ctg2+"\t"+dict_spe_edge_scaff[spe][adj].oriC1+"\t"+dict_spe_edge_scaff[spe][adj].oriC2+"\t"+str(dict_spe_edge_scaff[spe][adj].gap)+"\t"+dict_spe_edge_scaff[spe][adj].gf1+"\t"+dict_spe_edge_scaff[spe][adj].gf2+"\t"+dict_spe_edge_scaff[spe][adj].g1+"\t"+dict_spe_edge_scaff[spe][adj].g2+"\t"+dict_spe_edge_scaff[spe][adj].oriG1+"\t"+dict_spe_edge_scaff[spe][adj].oriG2+"\t"+str(dict_spe_edge_scaff[spe][adj].dist)+"\t"+str(dict_spe_edge_scaff[spe][adj].vscore)+"\t"+str(dict_spe_edge_scaff[spe][adj].dscore)+"\t"+str(dict_spe_edge_scaff[spe][adj].link)+"\n")
   scaff_gene_file.close()
   print "DONE"

   print "\n=> "+str(nb_scaff_kept)+"/"+str(Nb_scaff_edge_tot)+" edges have been kept (Others edges have at least one of the 2 linked contigs that have no gene in gene trees and have been erased OR linked several times the same couple fo genes (Keep the ADJ with best score))"

end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))

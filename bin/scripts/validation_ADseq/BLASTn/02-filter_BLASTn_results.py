#! /usr/bin/env python
# -*- coding: utf-8 -*-
### 
###   Goal:
###      Compare annotation gene files between original and after mapping of genes on genome assembly computed by minia
###         => RUN 6 experiments: 2 reads sampling (50% & ALL) x 3 species (Aalb, Aara & Adir)
###
###   INPUT:
###      1- INPUT Original GENE file
###         (data/data_DeCoSTAR/GENE_file)
###      2- BLATSn results file (Tab format)
###         (data/validation_ADseq/BLASTn/Anopheles_arabiensis/50pourc/align_Aara_ORI_50pourc_allparams.tab)
###      3- OUTPUT CTG file for MINIA genome assembly
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_arabiensis/50pourc/CTG_Aara_50pourc)
###      4- OUTPUT annotation file for MINIA genome assembly
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_arabiensis/50pourc/annotG_Aara_50pourc)
###      5- OUTPUT file for association between new and old CTG ID in MINIA genome assembly
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_arabiensis/50pourc/assoCTG_Aara_50pourc)
###      6- OUTPUT file containing list of genes unmapped on MINIA CTG
###         (data/validation_ADseq/ADseq/BLASTn/Anopheles_arabiensis/50pourc/unmapG_Aara_50pourc)
###      7- Species name
###         (Anopheles_arabiensis)
###      8- Minimum percentage identity
###         (90)
###      9- Minimum percentage coverage of query (CTG minia) on subject (CTG ori)
###         (90)
###      10- Verbose
###         (True/true/T/t or False/false/F/f)
###      11- Original minia assembly (FASTA file)
###         (data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_arabiensis/50pourc/minia_k59_m3_Aara_50pourc.contigs.fa)
###      12- OUTPUT New FASTA file with MERGED Minia contigs
###         (data/validation_ADseq/FASTA/SCAFF/blastn/Anopheles_arabiensis/minia_blastn_Aara_50pourc.fa)
###
###   OUTPUT: (RUN in ~1min15s)
###      - Create annotation gene file for MINIA genome assembly (param 3)
###      - Produce geneome assembly FASTA file to execute BESST and produce scaffolding adjacencies on BLASTn contigs/scaffolds 
###
###   Name: 02-filter_BLASTn_results.py           Author: Yoann Anselmetti
###   Creation date: 2016/07/27                   Last modification: 2017/02/27
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
from Bio import SeqIO

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


def ori_ctg(ori,seq_ctg,start,end):
   seq=""
   if ori=="-":
      # Reverse sequence
      seq_rev=rev(seq_ctg)
      size=len(seq_rev)
      # print size
      # Take start (that is close to the end of contig) and end (close to start of contig)
      start_rev=size-start
      end_rev=size-end
      seq=seq_rev[start_rev:end_rev+1]
   elif ori=="+":
      seq=seq_ctg[start-1:end]
   else:
      exit("ERROR, orientation should be \"+\" or \"-\" and not: "+ori)
   return seq

def rev(seq):
   return str(seq)[::-1]


def store_FASTA(FASTA_file):
   ctg_seq=dict()
   for ctg in SeqIO.parse(FASTA_file,"fasta"):
      ID=ctg.id.split(" ")[0]
      ctg_seq[ctg.id]=ctg.seq
   return ctg_seq


def write_files(select_spe,start_min,end_max,CTG_ID,size_max,CTG_list_on_cur_gene,GENE_list_on_cur_ctg,output_ctg,output_gene,output_assoctg,list_OLD_CTG_kept,list_NEW_CTG,ctg_seq_minia,seq_ctg_ref,output_fasta):
   ID=""
   seq_scaff=""
   end1=0
   ##########
   ### Create association file between OLD MINIA CTG and NEW MINIA CTG (merged)  AND  get the seq of the merged MINIA contig and write it in new FASTA file
   ##########
   minim=1000000000000000
   maxim=0
   size_tot=0
   for ctgM in CTG_list_on_cur_gene:
      ori=ctgM.ori
      seq_ctg_minia=ctg_seq_minia[ctgM.id]
      start=ctgM.ostart
      size=ctgM.size
      start_seqM=ctgM.mstart
      end_seqM=ctgM.mend

      list_OLD_CTG_kept.append(ctgM.id)
      if not ID:
         ID=ctgM.ori+ctgM.id
      else:
         ID+=":"+ctgM.ori+ctgM.id

      # size_align=0
      # if start_seqM<end_seqM:
      #    size_align=end_seqM-start_seqM+1
      # else:
      #    size_align=start_seqM-end_seqM+1

      # print str(size)+" - "+str(size_align)


      ###########
      ### Get sequence of merged contigs and write them in new FASTA file.
      ###########
      # If it's not the first contig of the scaffold, check if current contig overlap or not previous ctg
      if end1!=0:
         start2=ctgM.ostart
         # If contig OVERLAP then cut start of 2nd contig and concatenante sequence of the 2 contigs
         if end1>=start2:
            overlap=end1-start2+1
            seq=ori_ctg(ori,seq_ctg_minia,start_seqM,end_seqM)
            seq_scaff+=seq[overlap:]
            # TO DO?: Method to check if overlap is the same sequence for thez 2 contigs?

         # If contigs don't overlap then had X "N"'s nt in the gap between the 2 contigs
         else:
            gap=seq_ctg_ref[end1:start2-1]
            # print "\t"+ctgM.id
            # print "\tSEQ: "+str(end1+1)+" - "+str(start2-1)+"\t=)> "+gap
            seq_scaff+=gap
            seq=ori_ctg(ori,seq_ctg_minia,start_seqM,end_seqM)
            seq_scaff+=seq
      # if first contig of scaffolds just add ctg seq to scaffold seq
      else:
         seq=ori_ctg(ori,seq_ctg_minia,start_seqM,end_seqM)
         seq_scaff+=seq
      end1=ctgM.oend

      if end1>maxim:
         maxim=end1
      if start<minim:
         minim=start

      # print ctgM.id
      # print "seq_align_size: "+str(len(seq))
      # print seq



   ID_CTG=CTG_ID.split("__")[0]+"__len__"+str(end_max-start_min+1)
   output_assoctg.write(ID_CTG+"\t"+ID+"\n")

   # print str(start_min)+" - "+str(end_max)+"   VS   "+str(minim)+" - "+str(maxim)
   # print str(len(seq_scaff))+" == "+str(end_max-start_min+1)

   output_fasta.write(">"+ID_CTG+"\n")
   output_fasta.write(str(seq_scaff)+"\n")

   list_NEW_CTG.append(ID_CTG)
   if verbose:
      print "\t"+ID_CTG+" <- "+ID
   gene_start_gf=""
   gene_start_pos=0
   gene_start_ori=""
   gene_start_id=""
   gene_gf=""
   gene_start=0
   gene_end=0
   gene_ori=""
   gene_id=""
   for g in GENE_list_on_cur_ctg:
      gene_gf=g.gf
      gene_start=str(g.start-start_min+1)
      gene_end=str(g.end-start_min+1)
      gene_ori=g.ori
      gene_id=g.id
      if not gene_start_id:
         gene_start_gf=g.gf
         gene_start_pos=str(g.start-start_min+1)
         gene_start_ori=g.ori
         gene_start_id=g.id
      output_gene.write(select_spe+"\t"+ID_CTG+"\t"+gene_gf+"\t"+gene_id+"\t"+gene_ori+"\t"+gene_start+"\t"+gene_end+"\n")
   output_ctg.write(select_spe+"\t"+ID_CTG+"\t"+str(end_max-start_min+1)+"\t"+str(len(GENE_list_on_cur_ctg))+"\t"+gene_start_gf+"\t"+gene_start_id+"\t"+gene_start_ori+"\t"+str(gene_start_pos)+"\t"+gene_gf+"\t"+gene_id+"\t"+gene_ori+"\t"+gene_end+"\n")


def Nx(x,list_size,assembly_size):
   sum_size=0
   Nx=0
   for size in sorted(list_size,reverse=True):
      sum_size+=size
      # print size
      if sum_size>=float(assembly_size/(100.0/float(x))):
         Nx=size
         break
   return Nx

################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   ORI_GENE_file=open(argv[1],'r')
   blastn_file=open(argv[2],'r')
   output_CTG_file=argv[3]
   output_GENE_file=argv[4]
   output_assoCTG_file=argv[5]
   output_unmapG_file=argv[6]
   select_spe=argv[7]
   identity=float(argv[8])
   min_cov=int(argv[9])
   verbose=False
   if (argv[10]=="T" or argv[10]=="True" or argv[10]=="t" or argv[10]=="true"):
      verbose=True
   elif (argv[10]=="F" or argv[10]=="False" or argv[10]=="f" or argv[10]=="false"):
      verbose=False
   else:
      exit("ERROR, parameter 10 should be : T/True/t/true or F/False/f/false and not: "+argv[10]+" !")
   FASTA_file_minia=argv[11]
   output_FASTA_file=argv[12]
   FASTA_file_ref=argv[13]

   #######
   ### Store initial MINIA contigs in ctg_seq_minia dict() 
   #######

   print "Store initial MINIA contigs...",
   ctg_seq_minia=store_FASTA(FASTA_file_minia)
   print "DONE"


   print "Store contigs of reference genome...",
   ctg_seq_ref=store_FASTA(FASTA_file_ref)
   print "DONE"


   OUTPUT_DIR=path.dirname(path.realpath(output_CTG_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   OUTPUT_DIR=path.dirname(path.realpath(output_GENE_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   OUTPUT_DIR=path.dirname(path.realpath(output_assoCTG_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   OUTPUT_DIR=path.dirname(path.realpath(output_FASTA_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   GENE=namedtuple("GENE",["start","end","gf","id","ori"])
   CTGmap=namedtuple("CTGmap",["msize","octg","osize","ident","cov","len","mm","go","mstart","mend","ostart","oend","evalue","score"])
   CTGminia=namedtuple("CTGminia",["ostart","oend","mstart","mend","id","size","ori","score"])



######
### 1/ STORE LIST OF GENES / CTG ID PRESENT IN REFERENCE GENOME ASSEMBLY OF SELECTED_SPECIES
######
   print "Store list of genes corresponding to contig in reference genome assembly...",
   oriCTG_geneList=dict()
   for gene in ORI_GENE_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",gene)
      if r:
         name_species=r.group(1)
         contig_ID=r.group(2)
         GF_ID=r.group(3)
         gene_ID=r.group(4)
         orientation=r.group(5)
         start=r.group(6)
         end=r.group(7)
         exon_nb=r.group(8)
         exon_pos=r.group(9)

         if name_species==select_spe:
            if not contig_ID in oriCTG_geneList:
               oriCTG_geneList[contig_ID]=list()
            oriCTG_geneList[contig_ID].append(GENE(int(start),int(end),GF_ID,gene_ID,orientation))
      else:
           exit("ERROR line: "+gene+" is bad written")
   ORI_GENE_file.close()
   print "DONE"


######
### 2/ FILTER 1 of BLASTN CONSISTING TO KEEP ONLY ALIGNEMENT WITH COVERAGE>=90% AND IDENTITY>=90%
######
   print "Filter1 of BLASTn to kepp only alignement with coverage>=90% and identity>=90%...",
   miniaCTG_oriCTG=dict()
   ctg_present_after_filt1=set()

   for align in blastn_file:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$',align)
      if r:
         qname=r.group(1)
         qlen=int(r.group(2))
         sname=r.group(3)
         slen=int(r.group(4))
         pident=float(r.group(5))
         qcovhsp=int(r.group(6))
         length=int(r.group(7))
         mismatch=int(r.group(8))
         gapopen=int(r.group(9))
         qstart=int(r.group(10))
         qend=int(r.group(11))
         sstart=int(r.group(12))
         send=int(r.group(13))
         evalue=float(r.group(14))
         bitscore=float(r.group(15))
         
         # Keep only CTG minia with identity >= $(identity) %
         if (pident>=identity):
            # Keep only CTG minia with coverage on reference assembly >= $(min_cov) %
            if (qcovhsp>=min_cov):
               ctg_present_after_filt1.add(qname)
               ctg_align=CTGmap(qlen,sname,slen,pident,qcovhsp,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore)
               if not qname in miniaCTG_oriCTG:
                  miniaCTG_oriCTG[qname]=list()
               if not miniaCTG_oriCTG[qname]:
                  miniaCTG_oriCTG[qname].append(ctg_align)
               else:
                  ctg_align_str=miniaCTG_oriCTG[qname][0]
                  if (ctg_align_str.cov<=ctg_align.cov and ctg_align_str.ident<=ctg_align.ident and ctg_align_str.evalue>=ctg_align.evalue and ctg_align_str.score<=ctg_align.score):
                     miniaCTG_oriCTG[qname].append(ctg_align)
                     if (ctg_align_str.evalue>ctg_align.evalue and ctg_align_str.score<ctg_align.score):
                        miniaCTG_oriCTG[qname].remove(ctg_align_str)
   blastn_file.close()
   print "DONE"


   #######
   ### Browse list of MINIA ctg that have been mapped on reference genome assembly to compute MINIA assembly size in bp and contigs and N50 in bp 
   #######
   size_minia_after_filt1=0
   list_size_ctg_after_filt1=list()
   for ctg in ctg_present_after_filt1:
      # print ctg
      size_ctg=int(ctg.split("__")[2])
      size_minia_after_filt1+=size_ctg
      list_size_ctg_after_filt1.append(size_ctg)
   N50_filt1=Nx(50,list_size_ctg_after_filt1,size_minia_after_filt1)
   print "\t=> There are "+str(len(ctg_present_after_filt1))+" MINIA contigs still present after filter1, representing "+str(size_minia_after_filt1)+" bp (N50 (bp): "+str(N50_filt1)+")"
   list_size_ctg_after_filt1[:]=[]





######################
### 3/ KEEP ONLY MINIA CONTIGS WITH AN UNIQUE OPTIMAL COVERAGE AND IDENTITY ALIGNMENT ON REFERENCE GENOME
######################
   print "Filter2 of BLASTn consisting to keep only alignement with an unique optimal alignment score (in coverage and identity)...",
   oriCTG_miniaCTG=dict()
   list_ctg_remove_in_filter2=list()
   for ctg in miniaCTG_oriCTG:
      if verbose:
         print ctg
      if len(miniaCTG_oriCTG[ctg])==1:
         ctgmap=miniaCTG_oriCTG[ctg][0]
         if not ctgmap.octg in oriCTG_miniaCTG:
            oriCTG_miniaCTG[ctgmap.octg]=list()
         start=0
         end=0
         # Invert contig if orientation is inverted compare to original/reference contig
         orientation=""
         if (ctgmap.ostart<ctgmap.oend):
            mstart=ctgmap.mstart
            mend=ctgmap.mend
            ostart=ctgmap.ostart
            oend=ctgmap.oend
            orientation="+"
         else:
            mstart=ctgmap.mend
            mend=ctgmap.mstart
            ostart=ctgmap.oend
            oend=ctgmap.ostart
            orientation="-"

         ctgminia=CTGminia(ostart,oend,mstart,mend,ctg,ctgmap.msize,orientation,ctgmap.score)
         oriCTG_miniaCTG[ctgmap.octg].append(ctgminia)
      else:
         if verbose:
            print "\t"+ctg+" has "+str(len(miniaCTG_oriCTG[ctg]))+" alignments with maximum score on reference genome assembly of species "+select_spe
         for align in miniaCTG_oriCTG[ctg]:
            if verbose:
               print "\t"+"\t",
               print align
         list_ctg_remove_in_filter2.append(ctg)
   miniaCTG_oriCTG.clear()
   print "DONE"


   size_minia_rm_filter2=0
   for ctg in list_ctg_remove_in_filter2:
      size_minia_rm_filter2+=int(ctg.split("__")[2])

   size_minia_after_filt2=size_minia_after_filt1-size_minia_rm_filter2
   # Remove ctg from list of ctg after filter1 to compute    
   for ctg in list_ctg_remove_in_filter2:
      ctg_present_after_filt1.remove(ctg)

   list_size_ctg_after_filt2=list()
   for ctg in ctg_present_after_filt1:
      size_ctg=int(ctg.split("__")[2])
      list_size_ctg_after_filt2.append(size_ctg)
   ctg_present_after_filt1.clear()

   N50_filt2=Nx(50,list_size_ctg_after_filt2,size_minia_after_filt2)
   print "\t=> There are "+str(len(list_size_ctg_after_filt2))+" MINIA contigs still present after filter2, representing "+str(size_minia_after_filt2)+" bp (N50 (bp): "+str(N50_filt2)+")"






###############
### 4/ MERGING MINIA CONTIGS OVERLAPPING THE SAME GENE 
##############
   print "Filter3 - Merging of Minia contigs with gene overlapping several contigs...",
   output_ctg=open(output_CTG_file,'w')
   output_gene=open(output_GENE_file,'w')
   output_assoctg=open(output_assoCTG_file,'w')
   output_unmap=open(output_unmapG_file,'w')
   output_fasta=open(output_FASTA_file,'w')
   output_ctg.write("species\tctg\tctg_size\tctg_gene_nb\t5'_gene_family\t5'_gene\torientation_5'_gene\tstart_5'_gene\t3'_gene_family\t3'_gene\torientation_3'_gene\tend_3'_gene\n")
   output_gene.write("species\tctg\tgene_family\tID\torientation\tstart\tend\n")
   output_assoctg.write("NEW_CTG_ID\tOLD_CTG_ID\n")

   list_OLD_CTG_kept=list()
   list_NEW_CTG=list()
   for ctg in sorted(oriCTG_geneList):
      if verbose:
         print ctg,
      start_min=1000000000
      end_max=0
      CTG_ID=""
      size_max=0
      CTG_list_on_cur_gene=list()
      GENE_list_on_cur_ctg=list()
      gene_mapped=False

      if not ctg in oriCTG_miniaCTG:
         if verbose:
            print "\t=> None CTG from MINIA genome assembly have been mapped on CTG: "+ctg+" of REFERENCE/ORIGINAL genome assembly of species "+select_spe+"!!!"
         for gene in oriCTG_geneList[ctg]:
            output_unmap.write(gene.id+"\n")
            if verbose:
               print "\t=> None CTG from MINIA genome assembly have been mapped on this GENE !!!"
      else:
         # print "\n"+ctg+":"
         # Get seq of reference contig to fill gap between gene if necessary
         seq_ctg_ref=ctg_seq_ref[ctg]
         if verbose:
            print ":"
         # Browse list of genes in current contig of ORIGINAL/REFERENCE genome assembly of "select_spe"
         for gene in oriCTG_geneList[ctg]:
            start=gene.start
            end=gene.end
            # If last MINIA contig DOESN'T overlap previous and current gene in ORIGINAL/REFERENCE genome assembly
            if (end_max!=0 and start>end_max and gene_mapped):
               write_files(select_spe,start_min,end_max,CTG_ID,size_max,CTG_list_on_cur_gene,GENE_list_on_cur_ctg,output_ctg,output_gene,output_assoctg,list_OLD_CTG_kept,list_NEW_CTG,ctg_seq_minia,seq_ctg_ref,output_fasta)
               # Reinitialize parameters cause minia CTG from previous gene didn't overlap on current gene. 
               start_min=1000000000
               end_max=0
               CTG_ID=""
               size_max=0
               CTG_list_on_cur_gene=list()
               GENE_list_on_cur_ctg=list()
               gene_mapped=False
               # if verbose:
               #    print "\t\tstart_min: "+str(start_min)+" | end_max: "+str(end_max)

            # If gene start is lower than contig start AND If gene end is higher than contig end
            if start<start_min:
               start_min=start
            if end>end_max:
               end_max=end
            if verbose:
               print "\tGENE "+gene.id+": "+str(start)+" "+str(end)

            # Browse minia CTG mapped on current gene of ORIGINAL/REFERENCE genome assembly of "select_spe" to fuse them if several map on the gene
            for miniaCTG in sorted(oriCTG_miniaCTG[ctg]):
               if (miniaCTG.ostart<=end):
                  if (miniaCTG.oend>=start):
                     gene_mapped=True
                     if (miniaCTG.size>size_max):
                        CTG_ID=miniaCTG.id
                        size_max=miniaCTG.size
                     if (miniaCTG.ostart<start_min):
                        start_min=miniaCTG.ostart
                     if (miniaCTG.oend>end_max):
                        end_max=miniaCTG.oend
                     if not miniaCTG in CTG_list_on_cur_gene:
                        CTG_list_on_cur_gene.append(miniaCTG)
                     if verbose:
                        print "\t\tstart_min: "+str(start_min)+" | end_max: "+str(end_max)
               else:
                  break
            if not gene_mapped:
               output_unmap.write(gene.id+"\n")
               if verbose:
                  print "\t=> None CTG from MINIA genome assembly have been mapped on this GENE !!!"
            else:
               GENE_list_on_cur_ctg.append(gene)

         if gene_mapped:
            # When we arrive to the last gene of current contig print it.
            write_files(select_spe,start_min,end_max,CTG_ID,size_max,CTG_list_on_cur_gene,GENE_list_on_cur_ctg,output_ctg,output_gene,output_assoctg,list_OLD_CTG_kept,list_NEW_CTG,ctg_seq_minia,seq_ctg_ref,output_fasta)

   output_ctg.close()
   output_gene.close()
   output_assoctg.close()
   output_fasta.close()
   print "DONE"


   size_minia_after_filt3_OLD=0
   size_minia_after_filt3_NEW=0
   list_size_OLD_CTG_kept=list()
   list_size_NEW_CTG=list()
   for ctg in list_OLD_CTG_kept:
      size_ctg=int(ctg.split("__")[2])
      size_minia_after_filt3_OLD+=size_ctg
      list_size_OLD_CTG_kept.append(size_ctg)
   for ctg in list_NEW_CTG:
      size_ctg=int(ctg.split("__")[2])
      size_minia_after_filt3_NEW+=size_ctg
      list_size_NEW_CTG.append(size_ctg)
   list_OLD_CTG_kept[:]=[]
   list_NEW_CTG[:]=[]


   N50_filt3_OLD=Nx(50,list_size_OLD_CTG_kept,size_minia_after_filt3_OLD)
   N50_filt3_NEW=Nx(50,list_size_NEW_CTG,size_minia_after_filt3_NEW)

   print "\t=> There are "+str(len(list_size_OLD_CTG_kept))+" OLD MINIA contigs still present after filter3, representing "+str(size_minia_after_filt3_OLD)+" bp (N50 (bp): "+str(N50_filt3_OLD)+")"
   print "\t=> There are "+str(len(list_size_NEW_CTG))+" NEW MINIA contigs present after filter3, representing "+str(size_minia_after_filt3_NEW)+" bp (N50 (bp): "+str(N50_filt3_NEW)+")"
   list_size_OLD_CTG_kept[:]=[]
   list_size_NEW_CTG[:]=[]

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
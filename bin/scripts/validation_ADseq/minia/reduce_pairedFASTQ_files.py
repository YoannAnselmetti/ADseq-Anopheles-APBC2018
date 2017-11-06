#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from os import  path, makedirs, listdir, close, remove
from re import sub
from shutil import move
from tempfile import mkstemp
import random, itertools, errno, HTSeq, subprocess


 # fraction between 0 and 1
fraction = float(argv[1])
input_dir = argv[2]


#Method to change pattern in an other pattern in a file 
def replace(file_path, pattern, subst):
	"""
	Fonction to replace pattern in subst in file_path
	"""
	#Create temp file
	fh, abs_path = mkstemp()
	new_file = open(abs_path,'w')
	old_file = open(file_path)
	for line in old_file:
		new_file.write(line.replace(pattern, subst))
	#close temp file
	new_file.close()
	close(fh)
	old_file.close()
	#Remove original file
	remove(file_path)
	#Move new file
	move(abs_path, file_path)

def mkdir_p(dir_path):
	try:
		makedirs(dir_path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and path.isdir(dir_path):
			pass
		else:
			raise

def browse_dir(fraction,dir_path):
	if("SRR" in dir_path):
		SRR=path.basename(dir_path)	
		print "\nProcessing files: \n\t- "+dir_path+"/"+SRR+"_1.fastq.gz\n\t- "+dir_path+"/"+SRR+"_2.fastq.gz"
		reduce_FASTQ(fraction,dir_path+"/"+SRR+"_1.fastq",dir_path+"/"+SRR+"_2.fastq")
	else:
		for dir_file in sorted(listdir(dir_path)):
			browse_dir(fraction,dir_path+"/"+dir_file)

def reduce_FASTQ(fraction,file1,file2):

	out_file1 = sub('ALL','50pourc',file1)
	out_file2 = sub('ALL','50pourc',file2)

	OUTPUT_DIR=path.dirname(path.realpath(out_file1))
	# Create OUTPUT_DIR if not existing
	if not path.exists(OUTPUT_DIR):
		mkdir_p(OUTPUT_DIR)

	print "Creation of OUTPUT files:\n\t- "+out_file1+"\n\t- "+out_file2
	in1 = iter(HTSeq.FastqReader(file1+".gz"))
	in2 = iter(HTSeq.FastqReader(file2+".gz"))
	out1 = open(out_file1, "w")
	out2 = open(out_file2, "w")
	for read1, read2 in itertools.izip( in1, in2 ):
		if random.random() < fraction:
			read1.write_to_fastq_file(out1)
			read2.write_to_fastq_file(out2)
	out1.close()
	out2.close()

	print "\t\t=> gzip "+out_file1+"\n\t\t=> gzip "+out_file2
	command_line="gzip "+out_file1+"; gzip "+out_file2
	subprocess.call(command_line,shell=True)


browse_dir(fraction,input_dir)
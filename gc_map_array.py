#!/usr/bin/python
import numpy as np
import pandas as pd
import argparse
from six import iteritems
import sqlite3
import numpy.lib.recfunctions

# CREATE ARRAY WITH GC-BIAS AND MAPPABILITY INFORMATION
usage = '''This tool takes the human reference genome and a bedGraph-file with mappability scores (see http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability for 100 bp mappability in BigWig format, convert to bedGraph-format) as input and generates a csv file with gc-content and mappability score for every 100 bp''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--hg19', dest='hg19', help = 'Path to reference genome hg 19 directory', required= True)
parser.add_argument('--bed', dest='bedGraph', help = 'Path to BigWig file from http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability for 100 bp', required= True)
args = parser.parse_args()


def GC_MAP_ARRAY (hg19, bedGraph):
	path = hg19 + "/fasta/genome.fa"


	#map_dtype =  [('key', np.str), ('field', np.float), ('field', np.float),('field', np.float), ('field', np.float)]
	#gc_dtype = [('key', np.str), ('field', np.float)]
	gc_arr = np.empty((0,5), int)
	map_arr = np.empty((0,2), int)


	with open (path, 'r') as f_in:
		print 'Reads reference genome...'
		seq_dict = {}
		seq = []

		for line in f_in:
			if line[0] == ">":
				if len(seq) != 0: # if not first rount
					seq_string = ''.join(seq) # add seq as one string
					seq_dict[ID] = seq_string
					seq = [] # seq is now emty and ready for a new chromosome sequence

				ID_line = line[1:].split(' ')
				ID = ID_line[0].rstrip()
				print 'Processing: ', ID, '...'
				#print ID
			else:
				line = line.upper().strip() 	
				seq.append (line)

		# last chromosome		
		if len(seq) != 0: # if not first rount
			seq_string = ''.join(seq) # add seq as one string
			seq_dict[ID] = seq_string		

	window_size = 100
	
	with open (bedGraph, 'r') as bed_in:

		for key, value in seq_dict.iteritems(): 

			chrom = key
			print 'GC-content and mappability score for ', chrom, ' is prossesing...'
			seq_length = len(value)
			seq_list = []
			for nuc in value: # nake a list of all nucleotides in the genome for every chromosome specific
				seq_list.append(nuc)
			# count the GC content for every 100 bp
			for i in range (0, seq_length, window_size):

				pos_start = i+1
				pos_end = i+100
				seq_region = seq_list[i+1:i+100] 
				G = seq_region.count('G')
				C = seq_region.count('C') 
				GC_content = (G+C) / float(100)
				pos_start_digits = '%09d' % pos_start
				key_word = "%s:%s" % (chrom, pos_start_digits)

				gc_arr = np.append(gc_arr, np.array([[key_word, chrom, pos_start, pos_end, GC_content]]), axis=0)
		

			counter = 1
			for line in bed_in:
				line = line.split('\t')
				if line[0] == chrom:
					if int(line[1]) > counter:
						counter = int(line[2]) /100 * 100 +1
						#print counter			 
					
					region = True						
					while region == True:
						if counter >= int(line[1]) and counter <= int(line[2]):
							position_start = counter
							position_end = int(line[2])
							mappability_score = float(line[3])
							position_start_digits = '%09d' % position_start
							map_key_word = "%s:%s" % (chrom, position_start_digits)
							#print 'WE found the perfect macth', line[0], counter, counter + 99, line[3]
							map_arr = np.append(map_arr, np.array([[map_key_word, mappability_score]]), axis=0)
							counter += 100

						else:
							region = False	
							continue

	#print map_arr
	#print gc_arr 				
	data = gc_arr
	df1 = pd.DataFrame(data=data[1:,1:],
				index=data[1:,0], 
				#columns=data[0,1:])
				columns=['chrom', 'pos_start', 'pos_end', 'GC_content'])
	data2 = map_arr
	df2 = pd.DataFrame(data=data2[1:,1:],
				index=data2[1:,0], 
				columns=['mappability'])
	df = pd.concat([df1, df2], axis=1)
	conn = sqlite3.connect('normalization.db')

	df.to_sql("normalization_data", conn, if_exists="replace")
	#a = pd.read_sql_query("select * from normalization_data;", conn)
	#print a

array = GC_MAP_ARRAY(args.hg19, args.bedGraph)

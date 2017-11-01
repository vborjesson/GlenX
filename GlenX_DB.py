#!/usr/bin/python

import argparse
from six import iteritems
import sqlite3
import numpy as np
#import numpy.lib.recfunctions

#======================================================================================================================================
# GLENX.db
# Version 0.0.0 
#
# CREATE GLENX DATABASE 
# This is a script that creates a database containing information about GC-content and mappability score for the whole genome per 100 bp.
# It takes a bigGraph file generated from UCSCs BigWig-file, and calculates the mappability.score for every 100th bp in the genome. Then it takes 
# the human genome reference genome and calculats the GC-content(%) for every 100 bp. It is all saved in a SQLite database named GlenX.db
# 
# For more information: Visit https://github.com/vborjesson/GlenX.git
#======================================================================================================================================    

usage = '''This tool takes the human reference genome and a bedGraph-file with mappability scores (see http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability for 100 bp mappability in BigWig format, convert to bedGraph-format) as input and generates a csv file with gc-content and mappability score for every 100 bp''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--hg19', dest='hg19', help = 'Path to reference genome hg 19 directory', required= True)
parser.add_argument('--bed', dest='bedGraph', help = 'Path to BigWig file from http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability for 100 bp', required= True)
args = parser.parse_args()

def create_db (hg19, bedGraph):
	path = hg19 + "/fasta/genome.fa"
	valid_chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

	# Read in the reference genome file and calculate the GC-content per 100 bp
	with open (path, 'r') as f_in:
		print 'Reads reference genome...'
		seq_dict = {}
		seq = []
		not_valid = False

		for line in f_in:
			if line[0] == ">":

				if len(seq) != 0: # if not first rount
					seq_string = ''.join(seq) # add seq as one string
					seq_dict[ID] = seq_string
					seq = [] # seq is now emty and ready for a new chromosome sequence

				ID_line = line[1:].split(' ')
				ID = ID_line[0].rstrip() # chromosome
				if ID not in valid_chrom: # invalid chromosome	 -> continue to next!
					not_valid = True
					continue
				if ID in valid_chrom:
					not_valid = False
					print 'Processing: ', ID, '...'	
					continue					

			if not_valid == True:
				continue	
			else:
				line = line.upper().strip() 	
				seq.append (line)

		# last chromosome		
		if len(seq) != 0: # if not first rount
			seq_string = ''.join(seq) # add seq as one string
			seq_dict[ID] = seq_string		

	# Create a new database		
	db = sqlite3.connect('GlenX.db')
	cursorObject = db.cursor()
	cursorObject.execute('''CREATE TABLE GlenX (chrom TEXT NOT NULL, start_pos INT NOT NULL, end_pos INT, GC_content REAL, mappability_score REAL, PRIMARY KEY(chrom, start_pos));''')
	db.commit()

	for chrom, seq in seq_dict.iteritems():
		print 'GC-content for ', chrom, ' is prossesing...'
		seq_length = len(seq) / 100

		# make a list of all nucleotides in the genome for every chromosome specific
		seq_list = []
		for nuc in seq: 
			seq_list.append(nuc)

		# count the GC content for every 100 bp and add to SQlite database
		for x in range (0, 2500000):
			start_pos = x * 100

			end_pos = start_pos + 100
			if x > seq_length: # If there are no more information, just add 0 to database
				cursorObject.execute('''INSERT INTO GlenX (chrom, start_pos, end_pos, GC_content) VALUES (?,?,?,?)''', (chrom, start_pos, end_pos, 0))
				continue
			seq_region = seq_list[start_pos:end_pos] 
			G = seq_region.count('G')
			C = seq_region.count('C') 
			GC_content = (G+C) / float(100)		
			cursorObject.execute('''INSERT INTO GlenX (chrom, start_pos, end_pos, GC_content) VALUES (?,?,?,?)''', (chrom, start_pos, end_pos, GC_content))	
	db.commit()	

	with open (bedGraph, 'r') as bed_in:
		print 'Calculating mappability score for every 100 bp in the human reference genome...'
		score_arr = np.zeros((2500000, 100), float)
		prev_ch = ''
		new_chrom = False
		for line in bed_in:
			line = line.split('\t')
			if line[0] not in valid_chrom: # skip if it is not a valid chromosome
				continue
			ch = line[0] # chromosome
			if (ch != prev_ch) and (prev_ch != ''):  # new chromosome
				av_map_score = np.mean(score_arr, axis= 1)
				score_arr = np.zeros((2500000, 100), float) # Make new array for next chromosome
				for i in range (0, 2500000):
					start_pos = i * 100
					end_pos = start_pos + 100
					chrom = prev_ch
					mappability = av_map_score[i]
					cursorObject.execute('''UPDATE GlenX set mappability_score = ? WHERE chrom = ? AND start_pos = ? AND end_pos = ?''', (mappability, chrom, start_pos, end_pos))

			fr = int(line[1]) # start position with the same score
			to = int(line[2]) # end position with same score
			map_score = float(line[3])			
			for pos in range (fr, to):
				row = int (pos / 100)
				col = abs (pos) % 100
				score_arr[row, col] = map_score
			prev_ch = ch

		# add last chromosome to database
		av_map_score = np.mean(score_arr, axis= 1)
		for i in range (0, 2500000):
			start_pos = i * 100
			end_pos = start_pos + 100
			chrom = prev_ch
			mappability = av_map_score[i]
			cursorObject.execute('''UPDATE GlenX set mappability_score = ? WHERE chrom = ? AND start_pos = ? AND end_pos = ?''', (mappability, chrom, start_pos, end_pos))

		cursorObject.execute('''CREATE INDEX COV ON GlenX (chrom, start_pos, end_pos)''')	
		cursorObject.execute('''CREATE INDEX GC ON GlenX (GC_content)''')
		db.commit()
		db.close()
		print 'GlenX database completed'

array = create_db(args.hg19, args.bedGraph)

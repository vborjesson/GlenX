#!/usr/bin/python
import numpy as np
import argparse

# CREATE ARRAY WITH GC-BIAS AND MAPPABILITY INFORMATION
usage = '''Creates an file''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--hg19', dest='hg19', help = 'Path to reference genome hg 19 directory', required= True)
args = parser.parse_args()


def GC_MAP_ARRAY (hg19):
	path = hg19 + "/fasta/genome.fa"
	hg19_arr = np.empty((0,4), int)
	with open (path, 'r') as f_in:
		new_chrom = False
		counter = 0
		GC_count = 0
		for line in f_in:
			#line = line.upper()
			if line[0] == ">":
				new_chrom = True
				ID_line = line[1:].split(' ')
				ID = ID_line[0]
				#print ID
			else:
				line = line.upper().strip() 	
				for nuc in line:
					if new_chrom == True:
						counter = 0
						new_chrom = False
					counter += 1
					
        			if nuc == "G" or nuc == "C":
        				GC_count += 1
        			if counter % 100 == 0:
        				print ID, counter, GC_count

				'''
				if ID.startswith('NM_001289937.1'):
					print line.rstrip()
				'''
				window = 100

			#hg19_arr = np.append(hg19_arr, np.array([[ID, pos, GC, mappability]]), axis=0)	



array = GC_MAP_ARRAY(args.hg19)
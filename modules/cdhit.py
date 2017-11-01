#!/usr/bin/python
import os

# This module takes in an array containg sequenses that we want to cluster, creates a new fasta file and runs cdhit.
def ex_cdhit (array, id_path): # id_path = id_assembly/id_region_x_cdhit
	seq_id = 0
	fasta = '{}{}'.format(id_path, '.fasta')
	with open (fasta, 'w') as f_out:
		for line in array:
			seq_id += 1 
			sequence = line[10]
			f_out.write(seq_id)
			f_out.write(sequenses)
	cdhit_out = fasta.replace('.fasta', '')		
	process1 = ['cd-hit-v4.6.8-2017-0621/cd-hit-est -i', fasta, '-o', cdhit_out]
	os.system (' '.join(process1))


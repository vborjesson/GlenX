#!/usr/bin/python 

def tab_array (tab_file):
	import numpy as np
	from numpy import loadtxt, float32, dtype, int32
	import warnings	
	
	with open (tab_file, 'r') as tab_in:

		header = tab_in.readline()
		fields = header.split()

		fields_and_types = [('CHR','S16'), ('start',int), ('end',int),
                    ('coverage',float32), ('quality',float32)]
		data_dtype = dtype(fields_and_types)

		tab_array = loadtxt(tab_in, dtype=data_dtype,delimiter='\t', skiprows=1)

	np.seterr(divide='ignore', invalid='ignore')

	with warnings.catch_warnings():
		warnings.simplefilter("ignore", category=RuntimeWarning)
	
		chromosome_coverage = {}

		for a in range(1,23):
			a = str(a)
			chrom_num = tab_array['CHR'] == a
			cov_av = (tab_array['coverage'][chrom_num].mean())
			chromosome_coverage[a] = cov_av	

	return tab_array, chromosome_coverage	



	
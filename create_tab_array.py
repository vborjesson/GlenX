#!/usr/bin/python 

import numpy as np
from numpy import loadtxt, float32, dtype, int32

def tab_array (tab_file):	

	with open (tab_file, 'r') as tab_in:

		header = tab_in.readline()
		fields = header.split()

		fields_and_types = [('CHR','S16'), ('start',int), ('end',int), ('coverage',float32), ('quality',float32)]
		data_dtype = dtype(fields_and_types)

		tab_array = loadtxt(tab_in, dtype=data_dtype,delimiter='\t', skiprows=1)


	return tab_array	

#### SKIP THIS! 
'''
# get mean and standard deviation for read coverage over the entire chromosome
def normal_distr(tab_array, chromosome):
	chr_data = tab_array['CHR'] == chromosome
	cov_data = tab_array['coverage'][chr_data]
	data = np.delete(cov_data, 0)
	mu_chr_cov, std_chr_cov = norm.fit(data)
'''
'''
	mu_chr_cov = (tab_array['coverage'][data].mean())
	std_chr_cov = (tab_array['coverage'][data].std())
	
	return mu_chr_cov, std_chr_cov
'''	
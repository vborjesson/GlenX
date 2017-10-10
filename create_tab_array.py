#!/usr/bin/python 

import numpy as np
from numpy import loadtxt, float32, dtype, int32

def tab_array (tab_file):	

	with open (tab_file, 'r') as tab_in:
		print 'Creating array with read coverage information...'
		header = tab_in.readline()
		fields = header.split()

		fields_and_types = [('CHR','S16'), ('start',int), ('end',int), ('coverage',float32), ('quality',float32)]
		data_dtype = dtype(fields_and_types)

		tab_array = loadtxt(tab_in, dtype=data_dtype,delimiter='\t', skiprows=1)

		print 'Array completed'	
	return tab_array	

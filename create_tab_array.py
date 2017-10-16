#!/usr/bin/python 

import numpy as np
from numpy import loadtxt, float32, dtype, int32
import pandas as pd

def tab_array (tab_file):
	print 'Creataing a database containing read coverage information per 100 bp...'	
	#dtype_df =  {'CHR': np.str, 'start': np.int32, 'end': np.int32, 'coverage': np.float32, 'quality': np.float32, 'discordants': np.float32}
	tab_df = pd.read_csv(tab_file, sep='\t',  lineterminator='\n', names=['CHR', 'start', 'end', 'coverage', 'quality', 'discordants'], skiprows = 1, low_memory=False)#, usecols=['#CHR', 'start', 'end', 'coverage', 'quality'])
	tab_df['start'] = pd.to_numeric(tab_df['start'], errors='coerce')
	#tab_df['start'] = tab_df['start'].astype('int')
	#tab_df['coverage'] = tab_df['coverage'].astype('float32')

	'''
	with open (tab_file, 'r') as tab_in:
		print 'Creating array with read coverage information...'
		header = tab_in.readline()
		fields = header.split()

		fields_and_types = [('CHR','S16'), ('start',int), ('end',int), ('coverage',float32), ('quality',float32)]
		data_dtype = dtype(fields_and_types)

		tab_array = loadtxt(tab_in, dtype=data_dtype,delimiter='\t', skiprows=1)
	'''
	print 'dataframe completed'	
	return tab_df	

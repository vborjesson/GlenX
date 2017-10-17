#!/usr/bin/python 

import numpy as np
from numpy import loadtxt, float32, dtype, int32
import pandas as pd

def tab_array (tab_file):
	print 'Creataing a database containing read coverage information per 100 bp...'	

	tab_df = pd.read_csv(tab_file, sep='\t',  lineterminator='\n', names=['CHR', 'start', 'end', 'coverage', 'quality', 'discordants'], skiprows = 1, low_memory=False)#, usecols=['#CHR', 'start', 'end', 'coverage', 'quality'])
	tab_df['start'] = pd.to_numeric(tab_df['start'], errors='coerce')

	print 'dataframe completed'	
	return tab_df	

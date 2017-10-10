##!/usr/bin/python

import numpy as np
import pandas as pd
import sqlite3
import argparse
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

usage = '''Normalizign read coverage using mappability and GCcontent data from hg19 refernce genome''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--norm', dest='norm_db', help = 'Path to normalication database', required= True)
parser.add_argument('--tab', dest='tab', help = 'tab_array', required= True)
args = parser.parse_args()

def normalize_data(norm, tab):

	#conn = sqlite3.connect(norm)
	#df = pd.read_sql_query("SELECT * from normalization_data", conn)

	# normalize read depth from tab array using GC-bias and mappability. 
	data = tab 

	tab_df = pd.DataFrame(data=data[1:,1:],
				index=data[1:,0], 
				#columns=data[0,1:])
				columns=['chrom', 'start', 'end', 'read_coverage'])
	
	print (tab_df)
	tab_df['chrom'] = 'chr' + tab_df['chrom'].astype(str)
	tab_df['start'] = tab_df['start'].apply(lambda x: '{0:0>9}'.format(x))
	tab_df["key_name"] = tab_df["chrom"].map(str) + str(':') + tab_df["start"]
	
	#print (tab_df)
create = normalize_data(args.norm_db, tab)
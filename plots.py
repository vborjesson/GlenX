##!/usr/bin/python

import numpy as np
import pandas as pd
import sqlite3
import argparse
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

usage = '''Make plots for normalization data; mappability and GCcontent''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--norm', dest='norm_db', help = 'Path to reference genome hg 19 directory', required= True)
args = parser.parse_args()


conn = sqlite3.connect(args.norm_db)
df = pd.read_sql_query("SELECT * from normalization_data", conn)

gc_map_df = pd.DataFrame(df, columns=['GC_content', 'mappability'])
gc_map_df.plot.scatter(x='GC_content', y='mappability');


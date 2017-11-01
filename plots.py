##!/usr/bin/python

import numpy as np
import pandas as pd
import sqlite3
import argparse
import matplotlib.pyplot as plt
import matplotlib
#import ggplot
matplotlib.style.use('ggplot')

usage = '''Make plots for normalization data; mappability and GCcontent''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--norm', dest='norm_db', help = 'Path to reference genome hg 19 directory', required= True)
args = parser.parse_args()


conn = sqlite3.connect(args.db)

df = pd.read_sql_query("SELECT * from GlenX WHERE GC_content > 0 AND mappaility_score > 0", conn)

gc_map_df = pd.DataFrame(df, columns=['GC_content', 'mappability'])
p = gc_map_df.plot.scatter(x='GC_content', y='mappability');
ggsave(plot = p, filename='gc_mappability_plot.png')

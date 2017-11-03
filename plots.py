##!/usr/bin/python

import numpy as np
import pandas as pd
import sqlite3
import argparse
import matplotlib.pyplot as plt
import matplotlib
from ggplot import *

#import ggplot
matplotlib.style.use('ggplot')

usage = '''Make plots for normalization data; mappability and GCcontent''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--db', dest='db', help = 'Path to read_cov.db', required= True)
args = parser.parse_args()

db1 = sqlite3.connect(args.db)
db2 = sqlite3.connect('GlenX.db')
db1.execute('''ATTACH DATABASE 'GlenX.db' as db2''')
cursor = db1.cursor()

df = pd.read_sql_query("SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos and a.read_coverage between 1 and 100 and b.mappability_score > 0 LIMIT 20", db1)
print df

gc_map_df = pd.DataFrame(df, columns=['GC_content', 'read_coverage'])

p = gc_map_df.plot.scatter(x='GC_content', y='read_coverage');
ggsave(plot = p, filename='gc_ReadCoverage_plot.png')

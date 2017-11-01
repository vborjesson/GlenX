#/usr/bin/python

import sqlite3
import sys

db = sys.argv[1]
ch = 1
pos = 17900

def get_stats(db, ch, pos):

	statistics = {}
	# Connect the two databases
	db1 = sqlite3.connect(db)
	db2 = sqlite3.connect('GlenX.db')
	db1.execute('''ATTACH DATABASE 'GlenX.db' as db2''')
	cursor = db1.cursor()

	# get out the specific statistic.  
	gc_content = cursor.execute('''SELECT GC_content FROM db2.glen as a where a.chrom=? and a.start_pos = ? ''', (ch, pos)).fetchone()[0]
	r_i = cursor.execute('''SELECT read_coverage FROM read_cov where chrom=? and start_pos = ? ''', (ch, pos)).fetchone()[0]
	m_gc =cursor.execute('''SELECT avg(a.read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content FROM read_cov as a join db2.glen as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) a where a.GC_content=? and a.read_coverage > ?''', (gc_content,0)).fetchone()[0]
	m_all = cursor.execute('''SELECT avg(a.read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content FROM read_cov as a join db2.glen as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) a where a.read_coverage > ?''', (0,)).fetchone()[0]

	# normalize the read coverage for the breakpoint
	r_i_norm = r_i * m_all / m_gc

	statistics['r_i'] = r_i # read count for the i:th window
	statistics['m_all'] = m_all # average read count of all windows with a value higher then 0 
	statistics['m_gc'] = m_gc # average read count of all windows having the same gc percentage as i:th window
	statistics['gc_content'] = gc_content # gc content (refernece genome) for the i:th window
	statistics['r_i_norm'] = r_i_norm # normalized read counts for the i:th window

	db1.commit()

	return statistics

print statistics	



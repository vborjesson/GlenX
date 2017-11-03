#/usr/bin/python

import sqlite3
import sys

db = sys.argv[1]
ch = "1"
pos = 17900000

def get_stats(db, ch, pos):

	ch = '{}{}'.format('chr', ch)
	print 'Normalizing the read coverage and generating statistics for this breakpoint'

	statistics = {}
	# Connect the two databases
	db1 = sqlite3.connect(db)
	db2 = sqlite3.connect('GlenX.db')
	db1.execute('''ATTACH DATABASE 'GlenX.db' as db2''')
	cursor = db1.cursor()

	# get out the specific statistic.  
	map_i = cursor.execute('''SELECT mappability_score FROM db2.GlenX as a where a.chrom=? and a.start_pos = ? ''', (ch, pos)).fetchone()[0]
	gc_content = cursor.execute('''SELECT GC_content FROM db2.GlenX as a where a.chrom=? and a.start_pos = ? ''', (ch, pos)).fetchone()[0]
	r_i = cursor.execute('''SELECT read_coverage FROM read_cov where chrom=? and start_pos= ? ''', (ch, pos)).fetchone()[0]
	m_gc =cursor.execute('''SELECT avg(read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) as c where GC_content=? and read_coverage > ? and chrom=?''', (gc_content,0,ch)).fetchone()[0]
	m_all = cursor.execute('''SELECT avg(read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) as c where read_coverage > ? and chrom=?''', (0,ch)).fetchone()[0]
	m_all_at = cursor.execute('''SELECT avg(read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) as c where read_coverage > ? and chrom=? and mappability_score > ?''', (0,ch,0.9)).fetchone()[0]

	# normalize the read coverage for the breakpoint
	r_i_norm = r_i * m_all / m_gc

	statistics['r_i'] = r_i # read count for the i:th window
	statistics['m_all'] = m_all # average read count of all windows with a value higher then 0 
	statistics['m_gc'] = m_gc # average read count of all windows having the same gc percentage as i:th window
	statistics['gc_content'] = gc_content # gc content (refernece genome) for the i:th window
	statistics['r_i_norm'] = r_i_norm # normalized read counts for the i:th window
	statistics['map_i'] = map_i # mappability score for the ith window
	statistics['m_all_at'] = m_all_at # the average read coverage Above Thereshold
 
	db1.commit()

	return statistics

statistics = get_stats(db, ch, pos)
print statistics	
#print statistics['r_i']




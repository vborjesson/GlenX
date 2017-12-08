#/usr/bin/python

import sqlite3


#db = sys.argv[1]
#ch = "1"
#pos = 17900000

def get_stats(db, ch, pos):

	ch = ch.upper()
	ch = '{}{}'.format('chr', ch)
	print 'Normalizing the read coverage and generating statistics for this breakpoint'

	statistics = {}
	# Connect the two databases
	db1 = sqlite3.connect(db)
	db2 = sqlite3.connect('GlenX.db')
	db1.execute('''ATTACH DATABASE 'GlenX.db' as db2''')
	cursor = db1.cursor()

	print ch, pos
	
	# Classify SV using read depth over region
	if '|' in pos:
		pos = pos.split('|')
		pos1 = int(pos[0]) / 100 * 100
		pos2 = int(pos[1]) / 100 * 100

	# Classify SV using statistics for breakpoint
	else:	
		pos1 = int(pos) / 100 * 100
		pos2 = int(pos) + 100
		
	# get out the specific statistic.  
	cursor.execute('''SELECT mappability_score FROM db2.GlenX as a where a.chrom=? and a.start_pos = ? ''', (ch, pos1))
	map_i = cursor.fetchone()
	if map_i is None:
		return statistics
	else:
		map_i = map_i[0]	

	cursor.execute('''SELECT GC_content FROM db2.GlenX as a where a.chrom=? and a.start_pos between ? and ?''', (ch, pos1, pos2))
	gc_content = cursor.fetchone()
	if gc_content is None:
		return statistics
	else:
		gc_content = gc_content[0]	

	cursor.execute('''SELECT avg(read_coverage) FROM read_cov where chrom=? and start_pos between ? and ?''', (ch, pos1, pos2))
	r_i = cursor.fetchone()
	if r_i is None:
		return statistics
	else:
		r_i = r_i[0]	

	#else:	
	print gc_content
	m_gc =cursor.execute('''SELECT avg(read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) as c where GC_content=? and read_coverage > ? and chrom=?''', (gc_content,0,ch)).fetchone()[0]
	m_all = cursor.execute('''SELECT avg(read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) as c where read_coverage > ? and chrom=?''', (0,ch)).fetchone()[0]
	m_all_at = cursor.execute('''SELECT avg(read_coverage) FROM (SELECT a.chrom, a.start_pos, a.end_pos, a.read_coverage, b.GC_content, b.mappability_score FROM read_cov as a join db2.GlenX as b WHERE a.chrom = b.chrom and a.start_pos=b.start_pos) as c where read_coverage > ? and chrom=? and mappability_score > ?''', (0,ch,0.5)).fetchone()[0]

	# normalize the read coverage for the breakpoint
	print r_i, m_all_at, m_gc
	r_i_norm = r_i * m_all_at / m_gc

	statistics['r_i'] = r_i # read count for the i:th window
	statistics['m_all'] = m_all # average read count of all windows with a value higher then 0 
	statistics['m_gc'] = m_gc # average read count of all windows having the same gc percentage as i:th window
	statistics['gc_content'] = gc_content # gc content (refernece genome) for the i:th window
	statistics['r_i_norm'] = r_i_norm # normalized read counts for the i:th window
	statistics['map_i'] = map_i # mappability score for the ith window
	statistics['m_all_at'] = m_all_at # the average read coverage Above Thereshold

	db1.commit()

	return statistics






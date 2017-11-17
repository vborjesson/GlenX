#!/usr/bin/python 

import sqlite3
import os

#======================================================================
# Create a SQLite database containing read coverage information per 100 
# bp from tab file. 
#======================================================================	
def read_cov_db (ID, tab_file):
	print 'Creataing a database containing read coverage information per 100 bp...'	
	with open (tab_file, 'r') as f_in:
		db_name = '{}{}'.format(ID, '_read_cov.db')	
		if os.path.isfile(db_name) == False:
			db1 = sqlite3.connect(db_name)
			cursorObject = db1.cursor()
			cursorObject.execute('''CREATE TABLE read_cov (chrom TEXT NOT NULL, start_pos INT NOT NULL, end_pos INT, read_coverage REAL, PRIMARY KEY(chrom, start_pos));''')
			db1.commit()	
			for line in f_in:
				line = line.split('\t')
				chrom = '{}{}'.format('chr', line[0])
				start_pos = line[1]
				end_pos = line[2]
				cov = line[3]
				cursorObject.execute('''INSERT INTO read_cov (chrom, start_pos, end_pos, read_coverage) VALUES (?,?,?,?)''', (chrom, start_pos, end_pos, cov))
			cursorObject.execute('''CREATE INDEX REF ON read_cov (chrom, start_pos, end_pos)''')	
			db1.commit()
			db1.close()
	return db_name

#s=====================================================================
# Calculating and Returning median read coverage for the whole genome
#=====================================================================
def median_cov (db_path):
	db = sqlite3.connect(db_path)
	cursor_o = db.cursor()
	cursor_o.execute('''SELECT avg(read_coverage) FROM read_cov where read_coverage > ?''', (0,))
	med_cov = cursor_o.fetchone()[0]
	if med_cov is None:
		med_cov = 20
	db.commit()
	return med_cov
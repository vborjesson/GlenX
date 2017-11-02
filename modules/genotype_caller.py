##!/usr/bin/python

import numpy as np
import sys
import os
cwd = os.getcwd()
module_path = '{}{}'.format(cwd, '/modules')
sys.path.insert(0, module_path)
from check_bam_flag import bam_flag
from cigar import cigar_count
from cdhit import ex_cdhit
from statistics import get_stats


#===================================================================================================
# GENOTYPE CALLER
#===================================================================================================

def call_genotype (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, db):
	s_arr = np.empty((0,11), int) # Soft clipping alignment
	m_arr = np.empty((0,4), int) # matched alignment
	with open (sam, "r") as sam_in:

		#===========================================================================================
		#  Read in sam file and check for soft clips inside the specific region. 
		#===========================================================================================
		for line in sam_in:

			breakA = 0 # Breakpoint A (start) will be calculated below
			breakB = 0 # Breakpoint B (alt mapping) will be calculated below
			if line[0] == "@": # Skip info lines
				continue
			if line[0] == '\n': # if newline in end of file, skip this
				continue	

			else:
				line = line.upper().rstrip().split("\t")
				alt_chrA = str(line[2])
				if '.' in alt_chrA: # if the chromosome number contain a "." , it will be invalid and we will continue with next SV
					continue
				contig_start = int(line[3]) # start position for contig
				map_scoreA = int(line[4]) 
				cigar = line[5]
				strandA = bam_flag(line[1])
				
				if "S" in cigar:
					bad_quality = False # If there are several possible mate-mapping positions, this will be classified as bad quality and we will ignore these Breakpoints 
					SA = False # Second mapping position
					count_split_posA, cigar_length_posA = cigar_count (cigar, strandA)
					breakA += int(contig_start) # Breakpoint A  		
					breakA += count_split_posA 	

					# look at mate position of split reads. Can be found at optional field starting with SA:Z
					for field in line:
						if field.startswith("SA:"):
							split_info = field.split(":")
							positions = split_info[-1]
							n_position = positions.split(";") # split into number of positions. If more than one alternative position, skip! 
							if len(n_position) > 2: # due to one extra object; new line
								bad_quality = True
								break
							position = n_position[0].split(",")
							alt_chrB = str(position[0])
							mate_pos_start = position[1] 
							
							# strand 
							if position[2] == "+":
								strandB = 0
							elif position[2] == "-":
								strandB = 1
							map_scoreB = position[4]	
						
							count_split_posB, cigar_length_posB = cigar_count (position[3], strandB)
					
							breakB += int(mate_pos_start)
							breakB += count_split_posB
							SA = True							

						if field.startswith("AS:"):
							field = field.split(":")
							contig_l = field[-1]

					if bad_quality:
						continue
					if SA == False: # If the split contig have no second mapping place, continue
						continue			
					# count number of cigars, more cigars indicates untrustworthy SV. 
					cigar_length = 0
					cigar_length += cigar_length_posA
					cigar_length += cigar_length_posB	

					# check if breakpoints fall inside desired region
					region = False
					chromA = str(chromA)
					chromB = str(chromB)
					if alt_chrA == chromA and alt_chrB == chromB:
						if breakA >= posA_start and breakA <= posA_end and breakB >= posB_start and breakB <= posB_end:
							region = True	
						elif breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True	
					elif alt_chrA == chromB and alt_chrB == chromA:
						if breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True

					if region: # If region is True, save to array
						seq = line[8]
						s_arr = np.append(s_arr, np.array([[strandA, alt_chrA, breakA, map_scoreA, strandB, alt_chrB, breakB, map_scoreB, cigar_length, int(contig_l), seq]]), axis=0)	

				# Matched contig 		 
				elif "M" in cigar:
					count_match_pos, cigar_length_m = cigar_count (cigar, strandA) # count the number of base pairs that match to reference genome
					match_region_end = 0
					match_region_start = int(contig_start)
					match_region_end += match_region_start
					match_region_end += count_match_pos	
					m_arr = np.append(m_arr, np.array([[strandA, alt_chrA, match_region_start, match_region_end]]), axis=0)				

	# if no breakpoints could be found, skip this variant!					
	if len(s_arr) == 0:
		#print 'no breakpoints could be found'
		return s_arr, 'N/A', 'N/A' # returns N/A in order to continue the loop

	elif len(s_arr) > 0: # if one or more breakpoint was found		
	# If several predicted SVs in s_arr; compare them and choose the one that seems most trustworthy. 
	# Use map_score, cigar length, contig length and read coverage over breakpoints.
		if len(s_arr) > 1:
			# Best breakpoint will be the one with largest mapped contig. 
			seq_col = s_arr[:,10]
			seq_col = seq_col.astype(np.int)
			best_bp_number = np.argmax(seq_col)
			best_breakpoint = s_arr[best_bp_number]

		if len(s_arr) == 1:
			best_breakpoint = s_arr[0]		

		# Check if there is a matching contig to ref that will span over the predicted breakpoint. If there is; The genotype will be 
		# classified as heterozygous. If we can't find any matching contig spanning the breakpoint, we classify this as homozygous  
		genotype = ""
		for row in m_arr:
			if row[1] == best_breakpoint[1] and best_breakpoint[2] > row[2] and best_breakpoint[2] < row[3]:
				genotype = "0/1"
		if genotype == "":
			genotype = "1/1"

		# Get statistics from GlenX.db and read_cov.db
		statistics = get_stats(db, chromA, breakA)
		stat_map_score = statistics['map_i']

		# mappability threshold, we do not want to keep SVs who have a low mappability score = no support for SV.  
		#if stat_map_score is < 

		# classify SV. 	
		if best_breakpoint[1] != best_breakpoint[6]: # breakpoints are located on different chromosomes -> break end
			sv_type = "BND"
		if best_breakpoint[1] == best_breakpoint[6]: # breakpoints are located on the same chromosome
			if best_breakpoint[0] != best_breakpoint[5]: # sequences are in opposite directions -> inversed  
				sv_type = "INV"
			else: 
				sv_type = "improving code"	
	
		return best_breakpoint, genotype, sv_type, statistics


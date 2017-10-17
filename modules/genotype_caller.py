##!/usr/bin/python

import numpy as np

########################################### FUNCTION - GENOTYPE CALLER ########################################################################

def genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_df):
	s_arr = np.empty((0,13), int) # Soft clipping alignment
	m_arr = np.empty((0,4), int) # matched alignment
	with open (sam, "r") as sam_in:
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
				if '.' in alt_chrA: # if the chromosome number is longer than 2 letters, it will be invalid
					continue
				contig_start = int(line[3]) # start position for contig
				map_scoreA = int(line[4]) 
				cigar = line[5]
				strandA = bam_flag(line[1])
				
				if "S" in cigar:
					bad_quality = False # If there are several possible mate-mapping positions, this will be classified as bad quality and we will ignore these Breakpoints 
					SA = False# Second mapping position
					count_split_posA, cigar_length_posA = cigar_count (cigar, strandA)
					#print count_split_posA, cigar_length_posA
					breakA += int(contig_start) # Breakpoint A  		
					breakA += count_split_posA 	

					# look at mate position of split reads. Can be found at optional field starting with SA
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
						print 'bad quality, continuing with next SV'
						continue
					if SA == False: # If the split contig have no second mapping place, continue
						print 'No second mapping place have been predicted, continuing with next SV'
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
						print 'alt_chrA == chromA and alt_chrB == chromB is TRUE'
						if breakA >= posA_start and breakA <= posA_end and breakB >= posB_start and breakB <= posB_end:
							region = True	
						elif breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True	
					
					elif alt_chrA == chromB and alt_chrB == chromA:
						print 'alt_chrA == chromB and alt_chrB == chromA: is TRUE'
						if breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True

					if region:
						print 'region found'
						################################################
						################################################
						# average read depth for this position 100 bp
						# DO NOT FORGET TO CHANGE 1000 to 100
						##############################################
						##################################################

						# TIDDIT generates coverage information for every 100 bp. Therefore we will round down to nearest hundred and up to nearest hundred to get
						# specific region of interest, and calculate the average coverage for that region.  
						round_down_breakA = breakA/100
						round_down_breakA = int (round_down_breakA) * 100 # start position for region 100bp
						round_up_breakA = round_down_breakA + 100 # end position for region 100 bp

						cov1_breakA_list = tab_df.loc[(tab_df['CHR'] == str(alt_chrA)) * tab_df['start'] == str(round_down_breakA)]
						cov1_breakA = cov1_breakA_list['coverage'].astype(float)
						cov2_breakA_list = tab_df.loc[(tab_df['CHR'] == str(alt_chrA)) * tab_df['start'] == str(round_up_breakA)]
						cov2_breakA = cov2_breakA_list['coverage'].astype(float)
						av_cov_breakA = (cov1_breakA + cov2_breakA) / 2

						round_down_breakB = breakB/100
						round_down_breakB = int (round_down_breakB) * 100 # start position for region 100bp
						round_up_breakB = round_down_breakB + 100 # end position for region 100 bp
						
						cov1_breakB_list = tab_df.loc[(tab_df['CHR'] == str(alt_chrB)) * tab_df['start'] == str(round_down_breakB)]
						cov1_breakB = cov1_breakB_list['coverage'].astype(float)
						cov2_breakB_list = tab_df.loc[(tab_df['CHR'] == str(alt_chrB)) * tab_df['start'] == str(round_up_breakB)]
						cov2_breakB = cov2_breakB_list['coverage'].astype(float)
						av_cov_breakB = (cov1_breakB + cov2_breakB) / 2

						#print cov1_breakA, cov2_breakA, cov1_breakB, cov2_breakB
						#print 'cov break A, B', av_cov_breakA, av_cov_breakB	
						seq = line[8]
						s_arr = np.append(s_arr, np.array([[strandA, alt_chrA, breakA, map_scoreA, av_cov_breakA, strandB, alt_chrB, breakB, map_scoreB, av_cov_breakB, cigar_length, contig_l, seq]]), axis=0)	

						 
				elif "M" in cigar:
					count_match_pos, cigar_length_m = cigar_count (cigar, strandA) # count the number of base pairs that match to reference genome
					# print count_match_pos, cigar_length_m
					match_region_end = 0
					match_region_start = int(contig_start)
					match_region_end += match_region_start
					match_region_end += count_match_pos	
					m_arr = np.append(m_arr, np.array([[strandA, alt_chrA, match_region_start, match_region_end]]), axis=0)
	print s_arr
	print m_arr							

	# if no breakpoints could be found, skip this variant!					
	if len(s_arr) == 0:
		print 'no breakpoints could be found'
		return 'N/A', 'N/A', 'N/A'
		break

	if len(s_arr) > 0: # if one or more breakpoint was found		
	# If several predicted SVs in s_arr; compare them and choose the one that seems most trustworthy. 
	# Use map_score, cigar length, contig length and read coverage over breakpoints.
		if len(s_arr) > 1:
			Q_scoring = {}
			counter = 0
			for row in s_arr:
				# Standardize in range 0 - 60 and weight   
				Q_length = (((int(row[11]) - 0) * (60 - 0)) / (2000 - 0)) + 0
				Weight = 20
				Q_qscore = ((((int(row[3])+int(row[8])) - 0) * (60 - 0)) / (120 - 0)) + 0
				if int(row[10]) >= 5:
					if row[10] <= 7:
						Q_cigar = 60
					elif row[10] >= 8:
						Q_cigar = 120	
				else: 
					Q_cigar = 0
				Qscore = (Q_length * Weight) + Q_qscore - Q_cigar  		
				Q_scoring[counter] = Qscore		
				counter += 1
			best_breakpoint = max(Q_scoring, key=Q_scoring.get)		

		if len(s_arr) == 1:
			best_breakpoint = 0		


		# Check if there is a matching contig to ref that will span over the predicted breakpoint. If there is; The genotype will be 
		# classified as heterozygous. If we can't find any matching contig spanning the breakpoint, we classify this as homozygous  
		genotype = ""
		for row in m_arr:
			if row[1] == s_arr[best_breakpoint][1] and s_arr[best_breakpoint][2] > row[2] and s_arr[best_breakpoint][2] < row[3]:
				genotype = "0/1"
		if genotype == "":
			genotype = "1/1"

		sv_info = s_arr[best_breakpoint]
		#print sv_info

		# classify SV. 	
		if sv_info[1] != sv_info[6]: # breakpoints are located on different chromosomes -> break end
			sv_type = "BND"
		if sv_info[1] == sv_info[6]: # breakpoints are located on the same chromosome
			if sv_info[0] != sv_info[5]: # sequences are in opposite directions -> inversed  
				sv_type = "INV"
			else: 
				sv_type = "improving code"	

			#print sv_info[4], sv_info[9]   # compare average read coverage for that chromosome, and for breakpoint-region    	 	

		print 'genotype: ', genotype, 'type: ', sv_type	
		return sv_info, genotype, sv_type


#!/usr/bin/python

import argparse
#from GlenX import bam_flag
import math

usage = '''This function creates a new sam file containing only these SVs that are inside the region of interest''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--sam', dest='sam', help = 'path to x_splits.sam', required= True)
parser.add_argument('--out', dest='out', help = 'Path and name for outputfile (.sam)', required= True)
parser.add_argument('--region1', dest='region1', help = 'format chrom:start_pos-end_pos', required= True)
parser.add_argument('--region2', dest='region2', help = 'format chrom:start_pos-end_pos', required= True)

args = parser.parse_args()

sam = args.sam
out = args.sam
region1 = args.region1
region2 = args.region2

def check_region (sam, out, region1, region2):
	print 'function is initiated'
	with open (sam, 'r') as sam_in, open (out + '.sam', 'w') as f_out:
		for line in sam_in:
			print line

			breakA = 0 # Breakpoint A (start) will be calculated below
			breakB = 0 # Breakpoint B (alt mapping) will be calculated below
			if line[0] == "@": # Skip info lines
				continue
			if line[0] == '\n': # if newline in end of file, skip this
				continue	

			else:
				line = line.upper().rstrip().split("\t") # make a list of 
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

					# The specified regions we want our reads to fall within
					spec_regionA = region1.split(':').split('-')
					spec_regionB = region2.split(':').split('-')
					breakA = spec_regionA[0]
					posA_start = spec_regionA[1]
					posA_end = spec_regionA[2]
					breakB = spec_regionB[0]
					posB_start = spec_regionB[1]
					posB_end = spec_regionB[2]					

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

					string = ''.join(line)		
					if region:
						f_out.write(string)

print 'initiating check region..'
check_region (sam, out, region1, region2)


################ FUNCTION - BAM-FLAG converter ###############################

# Convert bam-flag to binary numbers and check if 2^4 (reverse strand) is present or not.
# Binary numbers with factor 2, backwards: 2^0, 2^1, 2^2, 2^3, 2^4 .. 2^n. In the list created below; if the 
# the first position (binary_list[0]) is 1, this means that 2^0 = 1 is present. So we will check if position 5
# (binary_list[4]) is 1 or 0, if it is 1 this means that 2^4 = 16 exist and strand is reverse.

def bam_flag (number):
	binary_list = [] 
	while number >= 1:
		number = float(number)
		number = number / 2
		#print number

		# If number is a whole number, this will be a binary 0. And if is a decimal number, it will be a 1. 
		if number.is_integer():
			binary_list.append(int(0))
		else:
			binary_list.append(int(1))
			# round number down to nearest integer
			number = math.floor(number)
	
	# Check if the fifth number (2^4) backwards in binary_list is a 1 
	if len(binary_list) >= 5:
		if binary_list[4] == 1:
			return 1
		else:
			return 0		
	else:
		return 0




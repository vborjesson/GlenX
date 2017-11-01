#!/usr/bin/python

import argparse
import sys 
import math
sys.path.insert(0, '/proj/b2014152/private/vanja/GlenX/modules')
from check_bam_flag import bam_flag
from cigar import cigar_count


usage = '''This function creates a new sam file containing only these SVs that are inside the region of interest''' 
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--sam', dest='sam', help = 'path to x_splits.sam', required= True)
parser.add_argument('--out', dest='out', help = 'Path and name for outputfile (.sam)', required= True)
parser.add_argument('--region1', dest='region1', help = 'format chrom:start_pos-end_pos', required= True)
parser.add_argument('--region2', dest='region2', help = 'format chrom:start_pos-end_pos', required= True)

args = parser.parse_args()

sam = args.sam
out = args.out
region1 = args.region1
region2 = args.region2

def check_region (sam, out, region1, region2):
	print 'function is initiated'
	with open (sam, 'r') as sam_in, open (out, 'w') as f_out:
		for line in sam_in:

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
						#print 'bad quality, continuing with next SV'
						continue
					if SA == False: # If the split contig have no second mapping place, continue
						#print 'No second mapping place have been predicted, continuing with next SV'
						continue			
					# count number of cigars, more cigars indicates untrustworthy SV. 
					cigar_length = 0
					cigar_length += cigar_length_posA
					cigar_length += cigar_length_posB	

					# check if breakpoints fall inside desired region
					region = False

					# The specified regions we want our reads to fall within

					spec_regionA = region1.split(':')
					spec_regionB = region2.split(':')
					chromA = spec_regionA[0]
					posA = spec_regionA[1]
					posA = posA.split('-')
					posA_start = int(posA[0])
					posA_end = int(posA[1])
					chromB = spec_regionB[0]
					posB = spec_regionB[1]
					posB = posB.split('-')
					posB_start = int(posB[0])
					posB_end = int(posB[1])					

					if alt_chrA == chromA and alt_chrB == chromB:
						if breakA >= posA_start and breakA <= posA_end and breakB >= posB_start and breakB <= posB_end:
							region = True	
						elif breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True	
					
					elif alt_chrA == chromB and alt_chrB == chromA:
						if breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True

					string = '\t'.join(line)			
					if region:
						print 'region is true'
						f_out.write('{} + {}' .format(string, '\n'))

		

print 'initiating check region..'
check_region (sam, out, region1, region2)







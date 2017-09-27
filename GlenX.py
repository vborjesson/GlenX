#!/usr/bin/python 

import sys
import os 
import argparse
import subprocess
import math
import numpy as np
from numpy import loadtxt, dtype,float32
import warnings
from create_tab_array import tab_array, normal_distr
from scipy.stats import norm

################# ARGPARSER ######################

usage = '''GlenX takes vcf- and bam-files as input and improve the prediction of genotype''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--vcf', dest='vcf_in', help = 'Path to vcf-file', required= False)
parser.add_argument('--bam', dest='bam_in', help = 'Path to bam-file', required= False)
parser.add_argument('--tab', dest='tab_in', help = 'Path to tab_file', required= False)
parser.add_argument('--ID', dest='ID', help= 'sample ID', required = False)
parser.add_argument('--sam', dest='sam', help= 'path to sam-file for dry run', required=False)
#parser.add_argument('--fa', dest= 'fa', help= 'Path to fasta-file with contigs generated from abyss', required = False)

args = parser.parse_args()

vcf = args.vcf_in
bam = args.bam_in
tab = args.tab_in
ID = args.ID
sam = args.sam
#fasta = args.fa

################# FUNCTION - REGION SPECIFIC ASSEMBLY (called variant-region) ##################################

# This function will generate a new bam file. For every called region; start-position +/- 1000bp and end-position +/- 1000bp are saved. 
# All reads in these regions are then extracted from the BAM-file into a new region-specific BAM-file. This function will call a bah script 
# that execute assembly, including converting bam to fasta, run assembly and mapping back to ref-genome. 
# This function will generate a sam-file. 


def region_specific_assembly (vcf, bam, ID, tab_array):
	
	counter = 0
	subprocess.call ('mkdir ' + ID + '_bam', shell = True)
	subprocess.call ('mkdir ' + ID + '_fasta', shell = True)
	subprocess.call ('mkdir ' + ID + '_assembly', shell = True)
	subprocess.call('chmod +x assembly.sh', shell=True)
	
	with open (vcf, "r") as vcf_in:
		for line in vcf_in:
			if line.startswith("#"):
				continue 
			else:
				split_line=line.lower().replace("chr","").split("\t")
				chromA = split_line[0]
				posA = split_line[1]
				posB = 0
				tags = split_line[7].split(";")
				for tag in tags:
					if tag.startswith("end="):
						tag = tag.split('=')
						chromB = chromA
						posB = tag[1]	
					elif tag.startswith("svtype=bnd"):
						alt_column = split_line[4].replace("n","").replace("[", "").replace("]","")
						alt_column = alt_column.split(":")
						chromB = alt_column[0]
						posB = alt_column[1]
				if posB == 0:
					continue	

				# unique ID for every region-specific bam-file	
				counter += 1	# Unique ID 
				region_ID = ID + "_region_" + str(counter)
				bam_file = ID + "_bam/" + region_ID + '.bam'
				fasta_file = ID + "_fasta/" + region_ID + '.fasta'
				sam_file = region_ID + '.sam'
				assembly_map = ID + '_assembly'

				posA = int(posA)
				posB = int(posB)
				posA_start = posA - 1000
				posA_end = posA + 1000
				posB_start = posB - 1000
				posB_end = posB + 1000	
				region2 = ""			  

				# Check for overlapping regions (we do not want double reads)
				if chromA == chromB:
					if posA < posB:
						if posA_end >= posB_start:
							region = str(chromA) + ':' + str(posA_start) + '-' + str(posB_end)
						else: 
							region = str(chromA) + ':' + str(posA_start) + '-' + str(posA_end) 
							region2 = str(chromB) + ':' + str(posB_start) + '-' + str(posB_end)

					if posA > posB:
						if posB_end >= posA_start:
							region = str(chromA) + ':' + str(posB_start) + '-' + str(posA_end)
						else:
							region = str(chromA) + ':' + str(posA_start) + '-' + str(posA_end) 
							region2 = str(chromB) + ':' + str(posB_start) + '-' + str(posB_end)				
				else:
					region = str(chromA) + ':' + str(posA_start) + '-' + str(posA_end)
					region2 = str(chromB) + ':' + str(posB_start) + '-' + str(posB_end)
				
				# average chromosome coverage divided in four, will give us a threshold of minimum read coverage that will be returned
				cov_med = np.median(tab_array['coverage'])
				print cov_med
				subprocess.call('./assembly.sh ' + bam + ' ' + region + ' ' + bam_file + ' ' + fasta_file + ' ' + sam_file + ' ' + assembly_map + ' ' + bwa_ref + ' ' + region_ID + ' ' + region2, shell = True)

	return chromA, chromB, posB_start, posA_end, posB_start, posB_end
	

########################################### FUNCTION - GENOTYPE CALLER ########################################################################

def genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_array):
	s_arr = np.empty((0,12), int) # Soft clipping alignment
	m_arr = np.empty((0,4), int) # matched alignment
	with open (sam, "r") as sam_in:
		for line in sam_in:
			breakA = 0 # Breakpoint A (start) will be calculated below
			breakB = 0 # Breakpoint B (alt mapping) will be calculated below
			if line[0] == "@": # Skip info lines
				continue
			else:
				line = line.upper().split("\t")
				alt_chrA = int(line[2])
				contig_start = int(line[3]) # start position for contig
				map_scoreA = int(line[4]) 
				cigar = line[5]
				strandA = bam_flag(line[1])
				
				if "S" in cigar:
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
							if len(n_position) > 2:
								break
							position = n_position[0].split(",")
							alt_chrB = int(position[0])
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

						if field.startswith("AS:"):
							field = field.split(":")
							contig_l = field[-1]

					# count number of cigars, more cigars indicates untrustworthy SV. 
					cigar_length = 0
					cigar_length += cigar_length_posA
					cigar_length += cigar_length_posB	

					# check if breakpoints fall inside desired region
					region = False
					if alt_chrA == chromA and alt_chrB == chromB:
						if breakA >= posA_start and breakA <= posA_end and breakB >= posB_start and breakB <= posB_end:
							region = True	
						elif breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True	
					
					elif alt_chrA == chromB and alt_chrB == chromA:
						if breakA >= posB_start and breakA <= posB_end and breakB >= posA_start and breakB <= posA_end:
							region = True

					if region:

						################################################
						################################################
						# average read depth for this position 100 bp
						# DO NOT FORGET TO CHANGE 1000 to 100
						##############################################
						##################################################

						# TIDDIT generates coverage information for every 100 bp. Therefore we will round down to nearest hundred and up to nearest hundred to get
						# specific region of interest, the calculate the average coverage for that region.  
						round_down_breakA = breakA/1000
						round_down_breakA = int (round_down_breakA) * 1000
						round_up_breakA = round_down_breakA + 1000
						chrom_breakA = chromA -1
						cov1_regionA = tab_arr['start'] == round_down_breakA
						cov1_breakA = (tab_arr['coverage'][cov1_regionA])	
						cov1_breakA = cov1_breakA[chrom_breakA]
						cov2_regionA = tab_arr['start'] == round_up_breakA
						cov2_breakA = (tab_arr['coverage'][cov2_regionA])	
						cov2_breakA = cov2_breakA[chrom_breakA]
						av_cov_breakA = (cov1_breakA + cov2_breakA) / 2

						round_down_breakB = breakB/1000
						round_down_breakB = int (round_down_breakB) * 1000
						round_up_breakB = round_down_breakB + 1000
						chrom_breakB = chromB -1
						cov1_regionB = tab_arr['start'] == round_down_breakB
						cov1_breakB = (tab_arr['coverage'][cov1_regionB])	
						cov1_breakB = cov1_breakB[chrom_breakB]
						cov2_regionB = tab_arr['start'] == round_up_breakB
						cov2_breakB = (tab_arr['coverage'][cov2_regionB])	
						cov2_breakB = cov2_breakB[chrom_breakB]
						av_cov_breakB = (cov1_breakB + cov2_breakB) / 2			

						#print cov1_breakA, cov2_breakA, cov1_breakB, cov2_breakB
						#print 'cov break A, B', av_cov_breakA, av_cov_breakB	
						
						s_arr = np.append(s_arr, np.array([[strandA, alt_chrA, breakA, map_scoreA, av_cov_breakA, strandB, alt_chrB, breakB, map_scoreB, av_cov_breakB, cigar_length, contig_l]]), axis=0)	
					
				elif "M" in cigar:
					count_match_pos, cigar_length_m = cigar_count (cigar, strandA) # count the number of base pairs that match to reference genome
					# print count_match_pos, cigar_length_m
					match_region_end = 0
					match_region_start = int(contig_start)
					match_region_end += match_region_start
					match_region_end += count_match_pos	
					m_arr = np.append(m_arr, np.array([[strandA, alt_chrA, match_region_start, match_region_end]]), axis=0)
	#print s_arr
	#print m_arr							

	# if no breakpoints could be found, skip this variant!					
	if len(s_arr) == 0:
		print 'no breakpoints could be found'
		return 'N/A'

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

	# Best scoring breakpoint prediction		
	best_breakpoint = max(Q_scoring, key=Q_scoring.get)

	# Check if there is a matching contig to ref that will span over the predicted breakpoint. If there is; The genotype will be 
	# classified as heterozygous. If we can't find any matching contig spanning the breakpoint, we classify this as homozygous  
	genotype = ""
	for row in m_arr:
		if row[1] == s_arr[best_breakpoint][1] and s_arr[best_breakpoint][2] > row[2] and s_arr[best_breakpoint][2] < row[3]:
			genotype = "0/1"
	if genotype == "":
		genotype = "1/1"

	sv_info = s_arr[best_breakpoint]
	print sv_info

	# classify SV. 	
	if sv_info[1] != sv_info[6]: # breakpoints are located on different chromosomes -> break end
		sv_type = "BND"
		print 'BND'
	if sv_info[1] == sv_info[6]: # breakpoints are located on the same chromosome
		mu_chr, std_chr = normal_distr(tab_arr, sv_info[1])
		if sv_info[0] != sv_info[5]: # sequences are in opposite directions -> inversed  
			sv_type = "INV"
		print mu_chr, std_chr	
		print sv_info[4], sv_info[9]   # compare average read coverage for that chromosome, and for breakpoint-region    	 	

	print genotype	

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


############################ FUNCTION - CIGAR counter #################################

# short function to count cigar field with soft clipping				
def cigar_count (cigar, strandA):
	n_bp = 0 # number of basepairs from start_pos
	for char in cigar:
		if char.isalpha():
			cigar = cigar.replace(char, char + ":")
			cigar_list = cigar.split(":")
			if strandA == 1: # reverse strand
				cigar_list = cigar_list[::-1] # flip list (complementary mapped to reference genome)	

	for cig in cigar_list:
		if "S" in cig:
			if n_bp == 0:
				cig = cig.replace("S", "")
				n_bp += int(cig)
				break
			else: 
				break  

		if "M" in cig:
			cig = cig.replace("M", "")
			n_bp += int(cig)

		if "D" in cig:
			cig = cig.replace("D", "")
			n_bp += int(cig)
	cigar_l = len(cigar_list) -1		
	return n_bp, cigar_l		



##### TEST DATA ######
chromA = 11
chromB = 11
posA_start = 26503276
posA_end = 26505276
posB_start = 30894624
posB_end = 30896624

# chromA, chromB, posB_start, posA_end, posB_start, posB_end = region_specific_assembly (vcf, bam, ID)

# Make an array from tab file containing read coverage information and average chromosome-specific coverage.  
tab_arr = tab_array(tab)	

#
# assembly = region_specific_assembly (vcf, bam, ID, tab_arr)

# Check if SVs are supported in de novo assemby, classify SV type and genotype.  
call_genotype = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_arr)


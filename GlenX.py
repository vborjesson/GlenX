#!/usr/bin/python 

import sys
import os 
import argparse
import subprocess
import math
import numpy as np

################# ARGPARSER ######################

usage = '''GlenX takes vcf- and bam-files as input and improve the prediction of genotype''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--vcf', dest='vcf_in', help = 'Path to vcf-file', required= False)
parser.add_argument('--bam', dest='bam_in', help = 'Path to bam-file', required= False)
parser.add_argument('--ID', dest='ID', help= 'sample ID', required = False)
parser.add_argument('--sam', dest='sam', help= 'path to sam-file for dry run', required=False)

args = parser.parse_args()

vcf = args.vcf_in
bam = args.bam_in
ID = args.ID
sam = args.sam

################# FUNCTION - REGION SPECIFIC ASSEMBLY (called variant-region) ##################################

# This function will generate a new bam file. For every called region; start-position +/- 1000bp and end-position +/- 1000bp are saved. 
# All reads in these regions are then extracted from the BAM-file into a new region-specific BAM-file. This function will call a bah script 
# that execute assembly, including converting bam to fasta, run assembly and mapping back to ref-genome. 


def region_specific_assembly (vcf, bam, ID):
	
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

				# uniqe ID for every region-specific bam-file	
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

				# Check for overlapping regions (we do not want doubble reads)
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
					
				subprocess.call('./assembly.sh ' + bam + ' ' + region + ' ' + bam_file + ' ' + fasta_file + ' ' + sam_file + ' ' + assembly_map + ' ' + bwa_ref + ' ' + region_ID + ' ' + region2, shell = True)


########################################### FUNCTION - GENOTYPE CALLER ########################################################################

def genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end):
	s_arr = np.empty((0,6), int) # Soft clipping alignment
	m_arr = np.empty((0,4), int) # matched alignment
	with open (sam, "r") as sam_in:
		for line in sam_in:
			if line[0] == "@":
				continue
			else:
				line = line.upper().split("\t")
				if "S" in line[5]:
					alt_chromA = int(line[2])
					alt_posA = int(line[3])
					alt_strandA = bam_flag(line[1])
					info = line[-1]
					info = info.split(":")
					alt = info[-1].split(",")
					alt_chromB = int(alt[0]) # alternative chromosome position 
					alt_posB = alt[1]
					alt_strandB = alt[2]
					if alt_strandB == '-':
						alt_strandB = 1
					elif alt_strandB == '+':
						alt_strandB = 0	 
					
					# print alt_chromA, chromA, alt_chromB, chromB
					
					# check if positions fall inside the desired region. If not; skip. 
					region = False

					if alt_chromA == chromA and alt_chromB == chromB:
						if alt_posA >= posA_start and alt_posA <= posA_end and alt_posB >= posB_start and alt_posB <= posB_end:
							region = True	
						elif int(alt_posA) >= posB_start and int(alt_posA) <= posB_end and int(alt_posB) >= posA_start and int(alt_posB) <= posA_end:
							region = True	
					
					elif alt_chromA == chromB and alt_chromB == chromA:
						if alt_posA >= posB_start and alt_posA <= posB_end and alt_posB >= posA_start and alt_posB <= posA_end:
							region = True	

					if region:
						# Look for breakpoints
						breakpointA = int(alt_posA)
						breakpointB = int(alt_posB)
						pre_breakpointA = line[5].split("S")
						pre_breakpointB = alt[3].split("S")
						bp_positionA = 0
						bp_positionB = 0

						# count how long in the contig the "S" (not match) is found. 
						# This is how it can look like 50S70M or 82M31S and the breakpoint is between M and S. 
						# One problem is that it also can be other alphabetic letters like 40M30S55H and so on. 

						for n in pre_breakpointA:
							if n.isdigit(): 
								if bp_positionA == 0: 
									breakpointA += int(n)
									break 
								else:
									breakpointA += bp_positionA
									break
							else:
								bp_positionA = ''.join([i for i in n if i.isdigit()])
								bp_positionA = int(bp_positionA)

						for v in pre_breakpointB:
							if v.isdigit():
								if bp_positionB == 0: 
									breakpointB += int(n)
									break 

								else:
									breakpointB += bp_positionB
							else:
								bp_positionB = ''.join([i for i in v if i.isdigit()])
								bp_positionB = int (bp_positionB)

						#print alt_strandA, breakpointA, alt_strandB, breakpointB
		
					#print line

						s_arr = np.append(s_arr, np.array([[alt_strandA, alt_chromA, breakpointA, alt_strandB, alt_chromB, breakpointB]]), axis=0)
				
				elif 'M' in line[5]:
					cigar = line[5].replace('M', '') # Check if only M cigar 
					if cigar.isdigit():
						match_contig = int(cigar)
						m_strand = bam_flag(line[1])
						m_chrom = int(line[2])
						m_pos_start = int(line[3])
						m_pos_end = 0
						m_pos_end += int(m_pos_start)
						m_pos_end += match_contig
						m_arr = np.append(m_arr, np.array([[m_strand,m_chrom,m_pos_start,m_pos_end]]), axis=0)
						#print match_contig, m_strand	

	# if no breakpoints count be found, skip this variant!					
	if len(s_arr) == 0:
		print 'no breakpoints could be found'
		return 'N/A' 

	#print m_arr, s_arr	

	# Check is there is a matching contig to ref, that will span over the predicted breakpoint. If there is; The genotype will be 
	# classified as heterozygot. If we cant find any matching contig spanning the breakpoint, we classify this as homozygous  
	genotype = ""
	for row in m_arr:
		if row[1] == s_arr[0][1] and s_arr[0][2] > row[2] and s_arr[0][2] < row[3]:
			genotype = "0/1"
	if genotype == "":
		genotype = "1/1"

	sv_info = s_arr[0]	

	# classify SV. 	
	if sv_info[1] != sv_info[4]:
		sv_type = "BND"
	if sv_info[1] == sv_info[4]:
		if sv_info[0] != sv_info[3]:
			sv_type = "INV"
		#elif 	 	

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
	# print binary_list[4]
	if len(binary_list) >= 5:
		if binary_list[4] == 1:
			return 1	
	else:
		return 0
				




chromA = 11
chromB = 11
posA_start = 26503276
posA_end = 26505276
posB_start = 30894624
posB_end = 30896624
#make_assembly = region_specific_assembly (vcf, bam, ID)
call_genotype = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end)



	
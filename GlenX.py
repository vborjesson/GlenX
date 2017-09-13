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
	arr = np.empty((0,3), int)
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
					
					# print alt_chromA, chromA, alt_chromB, chromB
					
					# check if positions fall inside the desired region. If not; skip. 
					region = False

					if alt_chromA == chromA and alt_chromB == chromB:
						if alt_posA >= posA_start and alt_posA <= posA_end and alt_posB >= posB_start and alt_posB <= posB_end:
							region = True	
						elif int(alt_posA) >= posB_start and int(alt_posA) <= posB_end and int(alt_posB) >= posA_start and int(alt_posB) <= posA_end:
							region = True	
					
					elif alt_chromA == chromB and alt_chromB == chromA:	
						print 'ja2'
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

						print alt_strandA, breakpointA, alt_strandB, breakpointB
		
					#print line

		arr = np.append(arr, np.array([[1,2,3]]), axis=0)
		arr = np.append(arr, np.array([[4,5,6]]), axis=0)
	print arr	



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

	if binary_list[4] == 1:
		return "-"
	if binary_list[4] == 0:
		return "+"	




chromA = 11
chromB = 11
posA_start = 26503276
posA_end = 26505276
posB_start = 30894624
posB_end = 30896624
#make_assembly = region_specific_assembly (vcf, bam, ID)
call_genotype = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end)



	
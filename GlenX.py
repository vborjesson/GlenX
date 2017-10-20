#!/usr/bin/python 

import sys
import os 
import argparse
import subprocess
import math
import numpy as np
from numpy import loadtxt, dtype,float32
import warnings
sys.path.insert(0, '/proj/b2014152/private/vanja/GlenX/modules')
from create_tab_array import tab_array
from genotype_caller import call_genotype
from cigar import cigar_count
from check_bam_flag import bam_flag
from scipy.stats import norm
import pandas as pd

################# ARGPARSER ######################

usage = '''GlenX takes vcf- and bam-files as input and improve the prediction of genotype''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--vcf', dest='vcf_in', help = 'Path to vcf-file', required= True)
parser.add_argument('--bam', dest='bam_in', help = 'Path to bam-file', required= True) 
parser.add_argument('--tab', dest='tab_in', help = 'Path to tab_file', required= True)
parser.add_argument('--norm_db', dest='norm_db', help = 'Path normalization.db', required= False)
parser.add_argument('--ID', dest='ID', help= 'sample ID', required = False)
parser.add_argument('--sam', dest='sam', help= 'path to sam-file for dry run', required=False)
parser.add_argument('--bwa_ref', dest='bwa_ref', help = 'Path to reference genome for bwa mem', required=True)
#parser.add_argument('--fa', dest= 'fa', help= 'Path to fasta-file with contigs generated from abyss', required = False)

args = parser.parse_args()
#sys.path.insert('./modules/')

vcf = args.vcf_in
bam = args.bam_in
tab = args.tab_in
ID = args.ID
sam = args.sam
bwa_ref = args.bwa_ref
#fasta = args.fa

################# FUNCTION - REGION SPECIFIC ASSEMBLY (called variant-region) ##################################

# This function will generate a new bam file. For every called region; start-position +/- 1000bp and end-position +/- 1000bp are saved. 
# All reads in these regions are then extracted from the BAM-file into a new region-specific BAM-file. This function will call a bah script 
# that execute assembly, including converting bam to fasta, run assembly and mapping back to ref-genome. 
# This function will generate a sam-file. 


def region_specific_assembly (vcf, bam, ID, tab_array, bwa_ref):	

	counter = 0
	subprocess.call ('mkdir ' + ID + '_bam', shell = True)
	subprocess.call ('mkdir ' + ID + '_fasta', shell = True)
	subprocess.call ('mkdir ' + ID + '_assembly', shell = True)
	subprocess.call ('mkdir ' + ID + '_GlenX_out', shell = True)
	subprocess.call('chmod +x assembly.sh', shell=True)

	file_name = '{}{}{}{}'.format(ID, '_GlenX_out/', ID, '_GlenX.vcf')
	
	with open (vcf, "r") as vcf_in, open (file_name, 'w') as f_out:
		print 'Reading VCF for region-specific assembly'
		ID_counter = 0
		for line in vcf_in:
			if line.startswith("#"):
				continue 
			else:
				print line
				split_line=line.lower().replace("chr","").split("\t")
				print split_line
				if len(split_line) <= 6: # This SV do not have all the correct information and will therefor not be analyzed 
					continue
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
				region_ID = '{}{}{}'.format(ID, "_region_", str(counter))
				assembly_map = '{}{}'.format(ID, '_assembly')
				

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
						if posA_end >= posB_start: # overlapping region
							region = '{}:{}-{}'.formart(str(chromA), str(posA_start), str(posB_end))
						else: 
							region = '{}:{}-{}'.format(str(chromA), str(posA_start), str(posA_end)) 
							region2 = '{}:{}-{}'.format(str(chromB), str(posB_start), str(posB_end))

					if posA > posB:
						if posB_end >= posA_start:
							region = '{}:{}-{}'.format(str(chromA), str(posB_start), str(posA_end))
						else:
							region = '{}:{}-{}'.format(str(chromA), str(posA_start), str(posA_end)) 
							region2 = '{}:{}-{}'.format(str(chromB), str(posB_start), str(posB_end))				
				else:
					region = '{}:{}-{}'.format(str(chromA), str(posA_start), str(posA_end))
					region2 = '{}:{}-{}'.format(str(chromB), str(posB_start), str(posB_end))
				
				# median chromosome coverage divided in four, will give us a threshold of minimum read coverage that will be returned
				cov_med = tab_arr['coverage'].median()# (tab_arr['coverage'])
				# minimum kmer coverage for de novo assembly
				min_cov = cov_med / 4
				print 'Minimum coverage accepted: ', min_cov
				print 'Initiate de novo assembly and mapping contigs back to reference'
				print 'region1-2', region, region2
				process = ['./assembly.sh', bam, region, ID, region_ID, str(min_cov), bwa_ref, region2]
				os.system(" ".join(process))
				sam = '{}/{}_mapped.sam' .format(assembly_map, region_ID)

			sv_info, genotype, sv_type = call_genotype (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_arr)	
			
			# Skipp this SV if no Breakpoint within the specified region could be find
			if sv_info == 'N/A':
				continue
			# Manipulate line with new improved SV information
			ID_counter += 1
			split_line[0] = chromA
			split_line[1] = sv_info[2] # Breakpoint for SV
			split_line[2] = 'SV_GlenX_{}'.format(str(ID_counter))
			split_line[4] = sv_type
			split_line[6] = 'PASS'
			if sv_type == 'BND':
			 	split_line[7] = 'SVTYPE=BND' #;CHRA=' + chromA  + ';CHRB=split_line[7] = 'SVTYPE=BND' #;CHRA=' + chromA  + ';CHRB=' + chromB + ';END=' + sv_info[7] ' + chromB + ';END=' + sv_info[7] # mating breakpoint 
				split_line[4] = 'N[{}:{}[' .format(chromB, sv_info[7])
			elif sv_type != 'BND':
				split_line[7] = 'SVTYPE={};CHRA={};CHRB={};END={}'.format(sv_type, chromA, chromB, sv_info[7]) #'SVTYPE=' + sv_type + ';CHRA=' + chromA  + ';CHRB=' + chromB + ';END=' + sv_info[7] 
			split_line[-1] = genotype	
			
			vcf_sv = '\t'.join(split_line) # make new tab seperated line of list, (preparations for writing to vcf-file) 
			string = str(vcf_sv)
			f_out.write(string)

	return 'KLART'		



##### TEST DATA ######
# chromA = 1
# chromB = 1
# posA_start = 26462
# posA_end = 28462
# posB_start = 30894
# posB_end = 32894

# chromA, chromB, posB_start, posA_end, posB_start, posB_end = region_specific_assembly (vcf, bam, ID)

# Make an array from tab file containing read coverage information and average chromosome-specific coverage.  
tab_arr = tab_array(tab)	
assembly = region_specific_assembly (vcf, bam, ID, tab_arr, bwa_ref)
#sv_info, genotype, sv_type = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_arr)	
# Check if SVs are supported in de novo assemby, classify SV type and genotype.  
# call_genotype = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_arr)


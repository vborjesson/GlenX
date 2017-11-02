#!/usr/bin/python 

import sys
import os 
import argparse
import subprocess
import math
import numpy as np
from numpy import loadtxt, dtype,float32
import warnings
cwd = os.getcwd()
module_path = '{}{}'.format(cwd, '/modules')
sys.path.insert(0, module_path)
from coverage_db import read_cov_db
from coverage_db import median_cov
from genotype_caller import call_genotype
from cigar import cigar_count
from check_bam_flag import bam_flag
from scipy.stats import norm
#import pandas as pd

################# ARGPARSER ######################

usage = '''GlenX takes vcf- and bam-files as input and improve the prediction of genotype''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--vcf', dest='vcf_in', help = 'Path to vcf-file', required= True)
parser.add_argument('--bam', dest='bam_in', help = 'Path to bam-file', required= True) 
parser.add_argument('--tab', dest='tab_in', help = 'Path to tab_file', required= True)
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


def region_specific_assembly (vcf, bam, ID, db, bwa_ref):	

	counter = 0
	subprocess.call ('mkdir ' + ID + '_bam', shell = True)
	subprocess.call ('mkdir ' + ID + '_fasta', shell = True)
	subprocess.call ('mkdir ' + ID + '_assembly', shell = True)
	subprocess.call ('mkdir ' + ID + '_GlenX_out', shell = True)
	subprocess.call('chmod +x assembly.sh', shell=True)

	file_name = '{}{}{}{}'.format(ID, '_GlenX_out/', ID, '_GlenX.vcf')
	
	with open (vcf, "r") as vcf_in, open (file_name, 'w') as f_out:
		print 'Reading VCF for region-specific assembly'
		ID_counter = 0 # all SVs will have I uniqe number
		info_field = False
		for line in vcf_in:
			if line.startswith("##"):
				f_out.write(line + "\n")
				if line.startswith("##INFO"):
					info_field = True 
				else: 
					if info_field == True:
						info_GlenX = "##INFO=<ID=GlenX,Number=6,Type=Float,Description='contig length, contig sequence, normalized breakpoint read-coverage (100bp), raw breakpoint read-coverage (100bp), gc-content(ref) and mappability-score(ref)'"
						info_field = False
				continue 
			else:
				split_line=line.lower().replace("chr","").split("\t")
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
				if posB == 0: # No secondary mapping have been found; continue with next SV
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
							region = '{}:{}-{}'.format(str(chromA), str(posA_start), str(posB_end))
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
				cov_med = median_cov (db)
				# minimum kmer coverage for de novo assembly
				min_cov = cov_med / 4
				print 'Minimum coverage accepted: ', min_cov
				print 'Initiate de novo assembly and mapping contigs back to reference'
				#print 'region1-2', region, region2
				process = ['./assembly.sh', bam, region, ID, region_ID, str(min_cov), bwa_ref, region2]
				os.system(" ".join(process))
				sam = '{}/{}_mapped.sam' .format(assembly_map, region_ID)

			sv_info, genotype, sv_type, statistics = call_genotype (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, db)	
			
			# If no breakpoints with good quality could be found inside our region, the old SV wiil be written to our new vcf
			if sv_type == 'N/A':
				f_out = write(line + '\n') 
				continue

			# Manipulate line with new improved SV information
			ID_counter += 1
			split_line[0] = chromA
			split_line[1] = sv_info[2] # Breakpoint for SV
			split_line[2] = 'SV_GlenX_{}'.format(str(ID_counter))
			split_line[4] = sv_type
			split_line[6] = 'PASS'
			
			contig_l = sv_info[9]
			contig_seq = sv_info[10]

			glenX_stats = '{}|{}|{}|{}|{}|{}'.format() # contig length, seq, normalized read ceverage, raw read-coverage, gc-content and mappabilty score
			
			if sv_type == 'BND':
			 	split_line[7] = 'SVTYPE=BND;GlenX={}'.format(GlenX_stats) #;CHRA=' + chromA  + ';CHRB=split_line[7] = 'SVTYPE=BND' #;CHRA=' + chromA  + ';CHRB=' + chromB + ';END=' + sv_info[7] ' + chromB + ';END=' + sv_info[7] # mating breakpoint 
				split_line[4] = 'N[{}:{}[' .format(chromB, sv_info[7])
			elif sv_type != 'BND':
				split_line[7] = 'SVTYPE={};CHRA={};CHRB={};END={}GlenX={}'.format(sv_type, chromA, chromB, sv_info[7], GlenX_stats) #'SVTYPE=' + sv_type + ';CHRA=' + chromA  + ';CHRB=' + chromB + ';END=' + sv_info[7] 
			split_line[-1] = genotype	
			
			vcf_sv = '\t'.join(split_line) # make new tab seperated line of list, (preparations for writing to vcf-file) 
			string = str(vcf_sv)
			f_out.write(string + "\n")

	return 'All variants in the vcf-file have been assembled, mapped and analyzed.'		


#======================================================================================================
# START GLENX
#======================================================================================================

db = read_cov_db(ID, tab)
print 'The database: ', db, 'is completed'	
assembly = region_specific_assembly (vcf, bam, ID, db, bwa_ref)
#sv_info, genotype, sv_type = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_arr)	
# Check if SVs are supported in de novo assemby, classify SV type and genotype.  
# call_genotype = genotype_caller (sam, chromA, chromB, posA_start, posA_end, posB_start, posB_end, tab_arr)


#!/usr/bin/python 

import sys
import os 
import argparse
import subprocess

################# ARGPARSER ######################

usage = '''GlenX takes vcf- and bam-files as input and improve the prediction of genotype''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--vcf', dest='vcf_in', help = 'Path to vcf-file', required= True)
parser.add_argument('--bam', dest='bam_in', help = 'Path to bam-file', required= True)
parser.add_argument('--ID', dest='ID', help= 'sample ID', required = True)

args = parser.parse_args()

vcf = args.vcf_in
bam = args.bam_in
ID = args.ID

################# FUNCTION - REGION SPECIFIC ASSEMBLY (called variant-region) ##################################

# This function will generate a new bam file. For every called region; start-position +/- 1000bp and end-position +/- 1000bp are saved. 
# All reads in these regions are then extracted from the BAM-file into a new region-specific BAM-file. This function will call a bah script 
# that execute assembly, including converting bam to fasta, run assembly and mapping back to ref-genome. 


def region_specific_assembly (vcf, bam, ID):
	
	counter = 0
	subprocess.call ('mkdir ' + ID + '_bam', shell = True)
	subprocess.call ('mkdir ' + ID + '_fasta', shell = True)
	subprocess.call('chmod +x assembly.sh', shell=True)
	
	with open (vcf, "r") as vcf_in:
		for line in vcf_in:
			if line.startswith("#"):
				continue 
			else:
				line=line.lower().replace("chr","").split("\t")
				chromA = line[0]
				posA = line[1]
				posB = 0
				tags = line[7].split(";")
				for tag in tags:
					if tag.startswith("end="):
						tag = tag.split('=')
						chromB = chromA
						posB = tag[1]	
					elif tag.startswith("svtype=bnd"):
						alt_column = line[4].replace("n","").replace("[", "").replace("]","")
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
					
				subprocess.call('./assembly.sh ' + bam + ' ' + region + ' ' + bam_file + ' ' + fasta_file + ' ' + ID + ' ' + region2, shell = True)

hej = region_specific_assembly (vcf, bam, ID)




	
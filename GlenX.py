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

################# FUNCTION - New BAM-file (called variant-region) ##################################

# This function will generate a new bam file. For every called region; start-position +/- 1000bp and end-position +/- 1000bp are saved. 
# All reads in these regions are then extracted from the BAM-file into a new region-specific BAM-file. 

def make_bam (vcf, bam, ID):
	counter = 0
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
				
				posA = int(posA)
				posB = int(posB)
				posA_start = posA - 1000
				posA_end = posA + 1000
				posB_start = posB - 1000
				posB_end = posB + 1000				  

				# Check for overlapping regions (we do not want doubble reads)
				if chromA == chromB:
					if posA < posB:
						if posA_end >= posB_start:
							subprocess.call('samtools view -b ' + bam + ' "Chr' + str(chromA) + ':' + str(posA_start) + '-' + str(posB_end) + '" > ' + ID + "_bam/" + region_ID + '.bam', shell = True)
							print 'a-b', posA, posB, chromA, posA_start, posB_end
							continue
					if posA > posB:
						if posB_end >= posA_start:
							subprocess.call('samtools view -b ' + bam + ' "Chr' + str(chromA) + ':' + str(posB_start) + '-' + str(posA_end) + '" > ' + ID + "_bam/" + region_ID + '.bam', shell= True)
							print "b-a", posA, posB, chromA, posB_start, posA_end
							continue					
				else:
					subprocess.call('samtools view -b ' + bam + ' "Chr' + str(chromA) + ':' + str(posA_start) + '-' + str(posA_end) + '"' + ' "Chr' + str(chromB) + ':' + str(posB_start) + '-' + str(posB_end) + '" > ' + ID + "_bam/" + region_ID + '.bam', shell=True) #, stdout=outfile)
					print posA, posB, chromA, posA_start, posA_end, chromB, posB_start, posB_end

subprocess.call ('mkdir ' + ID + '_bam', shell = True)
#subprocess.call ('samtools index ' + bam, shell = True)
hej = make_bam (vcf, bam, ID)




	
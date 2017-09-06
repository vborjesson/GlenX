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

args = parser.parse_args()

vcf = args.vcf_in
bam = args.bam_in

################# FUNCTION - New BAM-file (called variant-region) ##################################

# This function will generate a new bam file. For every called region; start-position +/- 1000bp and end-position +/- 1000bp are saved. 
# All reads in these regions are then extracted from the BAM-file into a new region-specific BAM-file. 

def region_bam (vcf, bam):
	with open (vcf, "r") as vcf_in:
		for line in vcf_in:
			if line.startswith("#"):
				continue 
			else:
				line=line.lower().replace("chr","").replace("end=","").split("\t")
				chrom = line[0]
				start = line[1]
				info_list = line[7].split(";")
				end = info_list[0]
				try:
   					val = int(end)
   					val = int(start)
				except ValueError:
   					print(start, end, "That's not an int!")
   					break
				start = int(start)
				end = int(end)
				start_a = start - 1000
				start_b = start + 1000
				end_a = end - 1000
				end_b = end + 1000

				print chrom, start_a, start_b, end_a, end_b 

				subprocess.call('samtools view -b '+ str(annotation_script), shell=True, stdout=outfile)
				print 'Annotation programs was added to the SVenX script'

hej = region_bam(vcf, bam)
	
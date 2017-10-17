#!/usr/bin/python


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
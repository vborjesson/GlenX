#!/usr/bin/python

################ FUNCTION - BAM-FLAG converter ###############################

# Convert bam-flag to binary numbers and check if 2^4 (reverse strand) is present or not.
# Binary numbers with factor 2, backwards: 2^0, 2^1, 2^2, 2^3, 2^4 .. 2^n. In the list created below; if the 
# the first position (binary_list[0]) is 1, this means that 2^0 = 1 is present. So we will check if position 5
# (binary_list[4]) is 1 or 0, if it is 1 this means that 2^4 = 16 exist and strand is reverse.

def bam_flag (number):
	import math
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
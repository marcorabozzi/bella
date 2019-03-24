#===================================================================================================
# Title:  Python script to create an identifier list from an overlapper output
# Author: G. Guidi
# Date:   14 Mar 2019
#===================================================================================================

# Overlappers' outputs need to be in PAF format: ``cd bench/'' and ``make paf''
# Run: ./id_list.py <overlapper-output-paf> <output-name>

import csv
import sys
import os

id_list = [] # overlapper output ids 	

with open(sys.argv[1], 'r') as file:
	# Here, we compare only read names, not alignment scores, etc.
	row = csv.reader(file, delimiter='\t')
	for row in file:
		row = row.rstrip().split("\t")
		id_list.append(row[0])	# paf format id1
		id_list.append(row[5])	# paf format id2

# The common approach to get a unique collection of items is to use a set. 
# Sets are unordered collections of distinct objects. 
# To create a set from any iterable, you can simply pass it to the built-in set() function.
# If you later need a real list again, you can similarly pass the set to the list() function.
# From: https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists.
id_list = list(set(id_list)) 		# removing duplicate entries 	

outfile = open(sys.argv[2], 'w')
outfile.write("\n".join(id_list))	# writing list to file	

print("End of program")

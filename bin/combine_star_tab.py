#!/usr/bin/env python3

import sys, csv

try:
	file1 = open(sys.argv[1])
except:
	sys.stderr.write('usage: script.py countsfile1 countsfile2 [...] > outfilename\n')
	sys.stderr.write('sums column 7 and 8 in all files based on matching the first 3 columns')
	sys.exit()

junction_info = {}
for fle in sys.argv[1:]:
	for line in open(fle):
		line = line.rstrip().split('\t')
		line[6] = int(line[6])  # unique counts
		line[7] = int(line[7])  # multimapper counts
		j = line[0]+':'+line[1]+'-'+line[2]
		if j not in junction_info:
			junction_info[j] = line
		else:
			junction_info[j][6] += line[6]
			junction_info[j][7] += line[7]

for j in junction_info:
	junction_info[j][6] = str(junction_info[j][6])
	junction_info[j][7] = str(junction_info[j][7])
	print('\t'.join(junction_info[j]))

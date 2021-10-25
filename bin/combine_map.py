#!/usr/bin/env python3

import sys, csv

s = 1
try:
	s = int(sys.argv[1])
	sys.argv = sys.argv[1:]
except:
	sys.stderr.write('usage: script.py n_supporting_reads mapfile1 mapfile2 [...] > outfilename\n')
	sys.stderr.write('combines isoform-read map files, remove entries with fewer than N supporting reads\n')
	sys.exit()

try:
	file1 = open(sys.argv[1])
except:
	sys.stderr.write('usage: script.py mapfile1 mapfile2 [...] > outfilename\n')
	sys.stderr.write('''combines isoform-read map files, could't open files\n''')
	sys.exit()

map_info = {}
for fle in sys.argv[1:]:
	for line in open(fle):
		iso, reads = line.rstrip().split('\t')
		if iso not in map_info:
			map_info[iso] = reads
		else:
			map_info[iso] += ','+reads

for i in map_info:
	if s == 1 or len(map_info[i].split(',')) >= s:
		print('\t'.join([i, map_info[i]]))

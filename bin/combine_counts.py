#!/usr/bin/env python3

import sys, csv, os

try:
	file1 = open(sys.argv[1])
except:
	sys.stderr.write('usage: combine_counts.py countsfile1 countsfile2 [...] outfilename\n')
	sys.exit()

counts = {}
salmon = False
for fle in sys.argv[1:-1]:
	for line in open(fle):
		line = line.rstrip().split('\t')
		if line[0] == 'Name':
			salmon = True
			continue
		if line[0] not in counts:
			counts[line[0]] = 0
		if salmon:
			counts[line[0]] += float(line[4])
		else:
			counts[line[0]] += float(line[1])

with open(sys.argv[-1], 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for transcript in counts:
		writer.writerow([transcript, counts[transcript]])
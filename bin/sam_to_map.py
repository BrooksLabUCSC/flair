#!/usr/bin/env python3
import sys, csv, os

try:
	sam = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py sam outfile\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for line in sam:
		line = line.rstrip().split('\t')
		if not line[0].startswith('@') and line[2] != '*':
			writer.writerow([line[0], line[2]])
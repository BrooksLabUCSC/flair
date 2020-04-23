#!/usr/bin/env python3
import sys, csv, os

try:
	psl = open(sys.argv[1])
	tsv = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py isoforms.psl counts_matrix.tsv isoforms_with_counts_output_name.psl\n')
	sys.exit()

tsv.readline()
counts = {}
for line in tsv:
	line = line.rstrip().split('\t')
	counts[line[0]] = line[1:]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for line in psl:
		line = line.rstrip().split('\t')
		if line[9] in counts:
			writer.writerow(line + counts[line[9]])
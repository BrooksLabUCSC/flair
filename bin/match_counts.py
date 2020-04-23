#!/usr/bin/env python3
import sys, csv, os

try:
	countsfile = open(sys.argv[1])
	psl = open(sys.argv[2])
	isbed = sys.argv[2][-3:].lower() != 'psl'
	min_reads = int(sys.argv[3])
	outfilename = sys.argv[4]  # sample_expression_pca_mtsf1.pdf
	append = len(sys.argv) > 5
except:
	sys.stderr.write('usage: script.py countsfile psl_to_append_counts min_read_threshold outfilename\n')
	sys.exit()

counts = {}
for line in countsfile:
	line = line.rstrip().split('\t')
	counts[line[0]] = float(line[1])

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for line in psl:
		line = line.rstrip().split('\t')
		name = line[3] if isbed else line[9]

		if name in counts:
			count = counts[name]
		else:
			count = 0
		if count >= min_reads:
			if append:
				writer.writerow(line + [count])
			else:
				writer.writerow(line)

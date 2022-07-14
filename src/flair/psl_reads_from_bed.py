#!/usr/bin/env python3
import sys, csv, os
keepnames = set()

for line in open(sys.argv[1]):  # bed of promoter-supported reads
	line = line.rstrip().split('\t')
	keepnames.add(line[3])

with open(sys.argv[3], 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	isbed = sys.argv[2][-3:].lower() != 'psl'
	for line in open(sys.argv[2]):  # psl or bed
		line = line.rstrip().split('\t')
		if isbed and line[3] in keepnames or not isbed and line[9] in keepnames:
			writer.writerow(line)

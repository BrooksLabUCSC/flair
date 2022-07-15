#!/usr/bin/env python3
import sys, csv, argparse, os

try:
	psl = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl'
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		nvrna = 'nvrna' in sys.argv[3]  # specify if stranded protocol
		reverse = 'reverse' in sys.argv[3]

	else:
		nvrna = reverse = False
except:
	sys.stderr.write('pull_starts.py psl|bed outfilename [nvrna]\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for line in psl:
		line = line.rstrip().split('\t')
		if isbed:
			chrom, name, start, end, strand, numblocks = line[0], line[3], int(line[1]), int(line[2]), line[5], line[9]
		else:
			chrom, name, start, end, strand, numblocks = line[13], line[9], int(line[15]), int(line[16]), line[8], line[17]

		if not nvrna and numblocks == 1:  # single exon gene, add both sides for cdna since strand is hard
			writer.writerow([chrom, start, start, name])
			writer.writerow([chrom, end, end, name])
			continue
		elif not reverse and '+' in strand or reverse and '-' in strand:
			tss = start
		elif not reverse and '-' in strand or reverse and '+' in strand:
			tss = end
		else:  # ambiguous strand, write both
			writer.writerow([chrom, start, start, name])
			tss = end	
		writer.writerow([chrom, tss, tss, name])
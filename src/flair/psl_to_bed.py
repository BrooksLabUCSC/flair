#!/usr/bin/env python3
import sys
import csv
import os
import bed

#try:
#	psl = open(sys.argv[1])
#	bed = sys.argv[2]
#except:
#	sys.stderr.write('usage: psl_to_bed in.psl out.bed\n')
#	sys.exit(1)
#
#with open(bed, 'wt') as outfile:
#	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
with open(sys.argv[1], 'r') as psl:
	for line in psl:
		line = line.rstrip().split('\t')
		if len(line) < 21:
			sys.stderr.write('fewer than 21 columns in the psl file, exiting\n')
			sys.exit(2)
		chrom, name, start, end = line[13], line[9], line[15], line[16]
		strand, blocksizes = line[8], line[18]
		bedObj = bed.Bed(chrom, start, end, name=name, strand=strand, numStdCols=12)
		bedObj.write(sys.stdin)
#		starts = line[20].split(',')[:-1]
#		relstarts = ','.join([str(int(n) - int(start)) for n in starts]) + ','
#		writer.writerow([chrom, start, end, name, '1000', strand, start, end, '0,0,0',
#			str(len(starts)), blocksizes, relstarts])

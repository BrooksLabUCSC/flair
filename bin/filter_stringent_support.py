#!/usr/bin/env python3
import sys, csv, os

try:
	isoforms = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl'
	alignment = open(sys.argv[2])
	minsupport = int(sys.argv[3])
	outfilename = sys.argv[4]
	if len(sys.argv) > 5:
		outfilename2 = sys.argv[5]
	else:
		outfilename2 = ''
	calculate_all = len(sys.argv) > 6
except:
	sys.stderr.write('usage: script.py isoforms.psl alignment.sam.psl minsupport out_isoforms.psl [out_assignments.txt] [calculate_all]\n')
	sys.exit(1)

isoform_info = {}
for line in isoforms:
	line = line.rstrip().split('\t')
	if isbed:
		blocksizes = [int(n) for n in line[10].split(',')[:-1]]
		name = line[3]
	else:
		blocksizes = [float(n) for n in line[18].split(',')[:-1]]
		name = line[9]
	isoform_info[name] = [sum(blocksizes), blocksizes[0], blocksizes[-1], line]

iso_read = {}  # isoform-read assignments for reads that span 25bp of the first and last exon
for line in alignment:  # reads aligned to the isoforms sam-turned-psl
	line = line.rstrip().split('\t')
	read, isoform = line[9], line[13]  # names

	if isoform not in iso_read:
		iso_read[isoform] = []
	elif len(iso_read[isoform]) > minsupport and not calculate_all:
		continue
	blocksizes = [int(n) for n in line[18].split(',')[:-1]]
	blockstarts = [int(n) for n in line[20].split(',')[:-1]]
	read_start, read_end = blockstarts[0], blockstarts[-1]+blocksizes[-1]

	info = isoform_info[isoform]
	isoform_length, first_blocksize, last_blocksize = info[0:3]

	right_coverage = left_coverage = False
	if len(blocksizes) == 1:  # single exon transcript
		if read_start < 25 and read_end > isoform_length - 25:
			right_coverage = left_coverage = True
	else:
		if first_blocksize < 25:
			if read_start < 2:
				left_coverage = True
		elif read_start <= (first_blocksize - 25):
			left_coverage = True

		if last_blocksize < 25:
			if (isoform_length - read_end) < 2:
				right_coverage = True
		if (isoform_length-last_blocksize + 25) <= read_end:
			right_coverage = True

	coverage = sum(blocksizes) / isoform_length
	# coverage = proportion of bases of the isoform that the read covers

	if right_coverage and left_coverage and coverage > 0.8:  
		iso_read[isoform] += [[read, isoform, coverage]]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for iso in iso_read:
		supporting = iso_read[iso]  # supporting reads
		if len(supporting) >= minsupport:
			writer.writerow(isoform_info[iso][3])

if outfilename2:  # map file
	with open(outfilename2, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
		for iso in iso_read:
			supporting = iso_read[iso]
			if len(supporting) >= minsupport:
				for s in supporting:
					writer.writerow(s)

#!/usr/bin/env python3

import sys, csv, os

try:
	counts_matrix = open(sys.argv[1])
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		sizefile = open(sys.argv[3])
	else:
		sizefile = ''
except:
	sys.stderr.write('usage: counts_to_tpm.py counts_matrix.tsv count_matrix.tpm.tsv [iso.sizes]\n')
	sys.stderr.write('convenience script for obtaining a file of isoform sizes: bin/fasta_seq_lengths.py\n')
	sys.stderr.write('if no isoform size file is provided, no length normalization will be done (just reads per million)\n')
	sys.exit()

sizes = {}
if sizefile:
	for line in sizefile:
		line = line.rstrip().split('\t')
		sizes[line[0]] = float(line[1])

header = counts_matrix.readline().rstrip().split('\t')
num_samples = len(header[1:])
matrix_data = [header]
all_rpk = [0] * num_samples
for line in counts_matrix:
	line = line.rstrip().split('\t')
	isoform_id = line[0]
	if sizes:
		rpk = [float(count)/sizes[isoform_id] for count in line[1:]]
	else:
		rpk = [float(count) for count in line[1:]]
	for n in range(num_samples):
		all_rpk[n] += rpk[n]
	matrix_data += [[isoform_id] + rpk]

all_rpk = [rpk/1e6 for rpk in all_rpk]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	writer.writerow(matrix_data[0])
	for line in matrix_data[1:]:
		for n in range(num_samples):
			line[n+1] = line[n+1]/all_rpk[n]
		writer.writerow(line)

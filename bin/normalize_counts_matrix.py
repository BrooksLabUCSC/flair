#!/usr/bin/env python3
import sys, csv, os
import numpy as np

try:
	matrix = open(sys.argv[1])
	outmatrix = sys.argv[2]
	if len(sys.argv) > 3:
		method = sys.argv[3].lower()
	else:
		method = 'cpm'
	if len(sys.argv) > 4:
		gtf = sys.argv[4]
	else:
		gtf = ''
except:
	sys.stderr.write('usage: script.py matrix outmatrix [cpm/uq/median] [gtf]\n')
	sys.stderr.write('gtf if normalization by protein coding gene counts only\n')
	sys.exit(1)

protein_coding = set()
if gtf:
	for line in open(gtf):
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
		if 'Y' in chrom or 'X' in chrom:
			continue
		if ty == 'gene':
			gene_type =  line[8][line[8].find('gene_type')+len('gene_type')+2:]
			gene_type = gene_type[:gene_type.find('"')]
			if gene_type == 'protein_coding':
				gene_id =  line[8][line[8].find('gene_id')+len('gene_id')+2:]
				gene_id = gene_id[:gene_id.find('"')]
				protein_coding.add(gene_id)

counts = {}  # header1: [0,1,1,0,2], header2: [0,3.6,0,2,1], ...
protein_coding_counts = {}
headers = matrix.readline().rstrip().split('\t')[1:]
for s in headers:
	counts[s] = []
	protein_coding_counts[s] = []

ids = []
for line in matrix:
	line = line.rstrip().split('\t')
	ids += [line[0]]
	for i in range(len(headers)):
		counts[headers[i]] += [float(line[i+1])]
	if protein_coding and line[0][line[0].rfind('_')+1:] in protein_coding:
		for i in range(len(headers)):
			protein_coding_counts[headers[i]] += [float(line[i+1])]

for sample in counts:

	if protein_coding:
		nonzero_counts = [c for c in protein_coding_counts[sample] if c != 0]
	else:
		nonzero_counts = [c for c in counts[sample] if c != 0]  # remove 0 counts
		
	if method == 'uq':
		scaling = np.percentile(nonzero_counts,75)
	elif method == 'median':
		scaling = np.median(nonzero_counts)
	else:
		scaling = sum(nonzero_counts)/1e6
	counts[sample] = [c / (scaling) for c in counts[sample]]  # normalize
	counts[sample] = [1 if c == 0 else c for c in counts[sample]]  # add pseudocount


with open(outmatrix, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	writer.writerow(['ids'] + headers)
	for i in range(len(ids)):  # i is line number
		row = [ids[i]]
		for j in range(len(headers)):  # j is header number
			row += [counts[headers[j]][i]]
		writer.writerow(row) 
	


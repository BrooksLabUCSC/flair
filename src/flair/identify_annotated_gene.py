#!/usr/bin/env python3
import sys, csv, os

try:
	psl = open(sys.argv[1])
	ref = open(sys.argv[2])
	outfilename = sys.argv[3]
	genepred = sys.argv[2][-3:].lower() == 'gp'
except:
	sys.stderr.write('usage: identify_annotated_gene.py psl ref.gtf/ref.gp isos_matched.psl \n')
	sys.exit(1)

def get_junctions(line):
	junctions = set()
	starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]  # for indexing purposes
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def bin_search(query, data):
	""" Query is a coordinate interval. Binary search for the query in sorted data, 
	which is a list of coordinates. Finishes when an overlapping value of query and 
	data exists and returns the index in data. """
	i = int(round(len(data)/2))  # binary search prep
	lower, upper = 0, len(data)
	while True:
		if upper - lower < 2:  # stop condition but not necessarily found
			break
		if data[i][1] < query[0]:
			lower = i
			i = int(round((i + upper)/2))
		elif data[i][0] > query[1]:
			upper = i
			i = int(round((lower + i)/2))
		else:  # found
			break
	return i

def contained(coords0, coords1, tol=0):
	""" complete coverage of coords0 by coords1, and coords0 can be tol larger.
	if coords0 is contained by coords1, then return the number of overlapping basepairs """
	if coords0[1] > coords1[0] and coords1[0] <= coords0[0]+tol and coords1[1] >= coords0[1]-tol:
		return min(coords1[1], coords0[1]) - max(coords1[0], coords0[0])
	return

prev_transcript, prev_exon = '', ''
all_juncs = {}  # matches a splice junction to gene name
all_se = {}  # single exon genes
# annotated_juncs = {}  # deprecated

if genepred:  # reading in annotated splice junctions
	for line in ref:
		line = line.rstrip().split('\t')
		gene, chrom, strand, numblocks = line[0], line[1], line[2], int(line[7])
		blockstarts = [int(n) + 1 for n in line[8].split(',')[:-1]]
		blockends = [int(n) for n in line[9].split(',')[:-1]]
		if chrom not in all_juncs:
			# annotated_juncs[chrom] = []
			all_juncs[chrom] = {}
			all_se[chrom] = []
		if numblocks == 1:
			all_se[chrom] += [(blockstarts[0], blockends[0])]
			continue
		# junctions = set()
		for start, end in zip(blockstarts[1:], blockends[:-1]):
			# junctions.add((end, start))
			all_juncs[chrom][(end, start)] = gene
		# annotated_juncs[chrom] += [(junctions, gene)]
else:
	for line in ref:  # extract all exons from the gtf, keep exons grouped by transcript
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
		if ty != 'exon':
			continue
		if chrom not in all_juncs:
			# annotated_juncs[chrom] = []
			all_juncs[chrom] = {}
			all_se[chrom] = []

		if 'gene_id' in line[8]:
			prev_gene = line[8][line[8].find('gene_id')+9:]
			prev_gene = prev_gene[:prev_gene.find('"')]
			this_transcript = line[8][line[8].find('transcript_id')+15:]
			this_transcript = this_transcript[:this_transcript.find('"')]
		elif 'geneid' in line[8].lower():
			prev_gene = line[8][line[8].find('geneid'):]
			prev_gene = prev_gene[:prev_gene.find(',')]
			this_transcript = line[8][line[8].find('transcript_id')+14:]
		else:
			sys.stderr.write('GTF format info column gene and transcript ids not recognized\n')
			sys.exit(1)

		if this_transcript != prev_transcript:
			if prev_transcript:
				if not junctions:  # single exon gene
					all_se[chrom] += [prev_exon]
			# 	annotated_juncs[chrom] += [(junctions, prev_gene)]
			junctions = set()
			prev_transcript = this_transcript
		elif strand == '-':
			junctions.add((end, prev_start))
			all_juncs[chrom][(end, prev_start)] = prev_gene
		else:
			junctions.add((prev_end, start))
			all_juncs[chrom][(prev_end, start)] = prev_gene
		prev_start = start
		prev_end = end
		prev_exon = (start, end, prev_gene)
	# annotated_juncs[chrom] += [(junctions, prev_transcript)]

for chrom in all_se:
	all_se[chrom] = sorted(list(all_se[chrom]), key=lambda x: x[0])

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for line in psl:
		line = line.rstrip().split('\t')

		if ';' in line[9][-3:]:
			line[9] = line[9][:line[9].rfind(';')]

		chrom = line[13]
		if chrom not in all_juncs:
			line[9] += '_chromnotinreference'
			writer.writerow(line)  # chrom not in the reference
			continue
		junctions = get_junctions(line)

		if '_EN' in line[9] or '_chr' in line[9] or '_chrom' in line[9]:  # already annotated
			writer.writerow(line)
			continue

		gene_hits = {}
		if not junctions:
			exon = (int(line[15]), int(line[16]))
			i = bin_search(exon, all_se[chrom])
			iscontained = False  # boolean: single-exon gene matches another s-e gene
			for e in all_se[chrom][i-2:i+2]:
				overlap = contained(exon, e, 20)
				if overlap:
					gene_hits[e[2]] = float(overlap)/(exon[1]-exon[0])  # gene name, % overlap
		else:
			for j in junctions:
				if j in all_juncs[chrom]:
					gene = all_juncs[chrom][j]
					if gene not in gene_hits:
						gene_hits[gene] = 0
					gene_hits[gene] += 1  # gene name, number of hits

		if not gene_hits:
			starts = [int(x) for x in line[20].split(',')[:-1]]  # block/exon starts
			line[9] += '_' + chrom + ':' + str(starts[0])[:-3] + '000'
			writer.writerow(line)
		else:
			genes = sorted(gene_hits.items(), key=lambda x: x[1])
			gene = genes[-1][0]
			line[9] += '_' + gene
			writer.writerow(line)


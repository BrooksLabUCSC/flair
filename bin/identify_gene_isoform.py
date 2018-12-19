import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
	
	isbed = sys.argv[1][-3:].lower() != 'psl' 
except:
	sys.stderr.write('usage: script.py psl/bed annotation.gtf renamed.psl/bed \n')
	sys.stderr.write('changes the name for each entry in psl/bed to the isoform and gene\n')
	sys.exit(1)

def get_junctions(line):
	junctions = set()
	starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]  # for indexing pupropses
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart + 1 for n in line[11].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[10].split(',')[:-1]]
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
junc_to_tn = {}  # chrom: {intron: [transcript_names], ... }
tn_to_juncs = {}  # chrom: {transcript_name: (junction1, junction2), ... }
all_se = {}
all_juncs = {}  # matches a splice junction to gene name

for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	if chrom not in junc_to_tn:
		junc_to_tn[chrom] = {}
		tn_to_juncs[chrom] = {}
		all_se[chrom] = []
		all_juncs[chrom] = {}
	this_transcript = line[8][line[8].find('transcript_id')+15:]
	this_transcript = this_transcript[:this_transcript.find('"')]

	if this_transcript != prev_transcript:
		if prev_transcript:
			if not junctions:  # single exon gene
				all_se[chrom] += [prev_exon]
			tn_to_juncs[chrom][prev_transcript] = junctions
			for j in junctions:
				if j not in junc_to_tn[chrom]:
					junc_to_tn[chrom][j] = set()
				junc_to_tn[chrom][j].add(prev_transcript)
		junctions = set()
		prev_transcript = this_transcript
	elif strand == '-':
		junctions.add((end, prev_start))
		all_juncs[chrom][(end, prev_start)] = prev_gene
	else:
		junctions.add((prev_end, start))
		all_juncs[chrom][(prev_end, start)] = prev_gene

	prev_start, prev_end = start, end
	prev_gene = line[8][line[8].find('gene_id')+9:]
	prev_gene = prev_gene[:prev_gene.find('"')]
	prev_exon = (start, end, prev_gene)

tn_to_juncs[chrom][this_transcript] = junctions
for j in junctions:
	junc_to_tn[chrom][j] = this_transcript

for chrom in all_se:
	all_se[chrom] = sorted(list(all_se[chrom]), key=lambda x: x[0])

novel, total = 0, 0
name_counts = {}  # to avoid redundant names
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		if isbed:
			junctions = get_junctions_bed12(line)
			chrom, name, start, end = line[0], line[3], int(line[1]), int(line[2])
		else:
			junctions = get_junctions(line)
			chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
		if chrom not in junc_to_tn:
			continue
		total += 1

		gene_hits = {}
		if not junctions:
			exon = (start, end)
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

		if not gene_hits:  # gene name will just be a chromosome locus
			gene = chrom + ':' + str(start)[:-3] + '000'
		else:  # gene name will be whichever gene the entry has more shared junctions with
			genes = sorted(gene_hits.items(), key=lambda x: x[1])
			gene = genes[-1][0]

		transcript = ''
		if junctions:
			matches = set()
			for j in junctions:
				if j in junc_to_tn[chrom]:
					matches.update(junc_to_tn[chrom][j])
			for t in matches:
				if tn_to_juncs[chrom][t] == junctions:
					transcript = t  # annotated transcript identified
					break

		if not transcript:
			novel += 1
			newname = name + '_' + gene
		else:
			newname = transcript + '_' + gene

		if newname not in name_counts:
			name_counts[newname] = 0
		else:
			name_counts[newname] += 1
			newname = newname + '-' + str(name_counts[newname])

		if isbed:
			line[3] = newname
		else:
			line[9] = newname
		writer.writerow(line)

# sys.stderr.write('{} out of {} isoforms have novel splice junction chains\n'.format(novel, total))


import sys, csv

try:
	psl = open(sys.argv[1])
	ref = open(sys.argv[2])
	outfilename = sys.argv[3]
	genepred = sys.argv[2][-3:] != 'gtf'
except:
	sys.stderr.write('usage: script.py psl ref.gtf/ref.gp isos_matched.psl \n')
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
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

prev_transcript = ''
annotated_juncs = {}

all_juncs = {}

if genepred:  # reading in annotated splice junctions
	for line in ref:
		line = line.rstrip().split('\t')
		gene, chrom, strand = line[0], line[1], line[2]
		blockstarts = [int(n) + 1 for n in line[8].split(',')[:-1]]
		blockends = [int(n) for n in line[9].split(',')[:-1]]

		if chrom not in annotated_juncs:
			annotated_juncs[chrom] = []
			all_juncs[chrom] = {}
		junctions = set()
		for start, end in zip(blockstarts[1:], blockends[:-1]):
			junctions.add((end, start))
			all_juncs[chrom][(end, start)] = gene
		annotated_juncs[chrom] += [(junctions, gene)]
else:
	for line in ref:  # extract all exons from the gtf, keep exons grouped by transcript
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
		if ty != 'exon':
			continue
		if chrom not in annotated_juncs:
			annotated_juncs[chrom] = []
			all_juncs[chrom] = {}
		this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
		this_transcript = line[8][line[8].find('transcript_id')+15:]
		this_transcript = this_transcript[:this_transcript.find('"')]

		prev_gene = line[8][line[8].find('gene_id')+9:]
		prev_gene = prev_gene[:prev_gene.find('"')]

		if this_transcript != prev_transcript:
			if prev_transcript:
				annotated_juncs[chrom] += [(junctions, prev_gene)]
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
	annotated_juncs[chrom] += [(junctions, prev_transcript)]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		chrom = line[13]
		if chrom not in annotated_juncs:
			continue
		junctions = get_junctions(line)

		if '_EN' in line[9]:  # already annotated
			writer.writerow(line)
			continue

		gene_hits = {}
		for j in junctions:
			if j in all_juncs[chrom]:
				gene = all_juncs[chrom][j]
				if gene not in gene_hits:
					gene_hits[gene] = 0
				gene_hits[gene] += 1

		if not gene_hits:
			starts = [int(x) for x in line[20].split(',')[:-1]]  # block/exon starts
			line[9] += '_' + chrom + ':' + str(starts[0])[:-3] + '000'
			writer.writerow(line)
		else:
			genes = sorted(gene_hits.items(), key=lambda x: x[1])
			gene = genes[-1][0]
			line[9] += '_' + gene
			writer.writerow(line)


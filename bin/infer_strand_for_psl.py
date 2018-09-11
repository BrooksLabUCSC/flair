import sys, csv

try:
	psl = open(sys.argv[1])
	ref = open(sys.argv[2])
	outfilename = sys.argv[3]
	genepred = len(sys.argv) > 4
except:
	sys.stderr.write('usage: script.py psl ref.gtf out.psl\n')
	sys.stderr.write('usage: script.py psl ref.gp out.psl \'gp\'\n')
	sys.exit()

annotmin, annotpos = {}, {}
if genepred:  # reading in annotated splice junctions
	for line in ref:
		line = line.rstrip().split('\t')
		gene, chrom, strand = line[0], line[1], line[2]
		blockstarts = [int(n) for n in line[8].split(',')[:-1]]
		blockends = [int(n) for n in line[9].split(',')[:-1]]
		for start, end in zip(blockstarts, blockends):
			if strand == '+':
				if chrom not in annotpos:
					annotpos[chrom] = {}
				annotpos[chrom][start] = gene
				annotpos[chrom][end] = gene
			else:
				if chrom not in annotmin:
					annotmin[chrom] = {}
				annotmin[chrom][start] = gene
				annotmin[chrom][end] = gene
else:
	for line in ref:
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8]
		if ty != 'exon':
			continue
		gene = gene[gene.find('gene_id')+len('gene_id')+2:]
		gene = gene[:gene.find('"')]
		if strand == '+':
			if chrom not in annotpos:
				annotpos[chrom] = {}
			annotpos[chrom][start] = gene
			annotpos[chrom][end] = gene
		else:
			if chrom not in annotmin:
				annotmin[chrom] = {}
			annotmin[chrom][start] = gene
			annotmin[chrom][end] = gene

def find_wiggle(coord, annot, annot2={}, maxdist=100):
	""" Finds the distance between coordinate and the closest annotated pos in annot dict. """
	wiggle = 0
	while coord + wiggle not in annot and coord + wiggle not in annot2:
		if wiggle == maxdist:
			break
		if wiggle == 0:
			wiggle += 1
		elif wiggle >= 0:
			wiggle = wiggle * -1
		else:
			wiggle = (wiggle-1) * -1
	return wiggle

def find_strand(annotpos, annotmin, jstarts, jends, maxdist, skipsize, skipnum=3):
	strands = []  # tally of predicted strands for each junction
	gene_candidates = []
	for start, end in zip(jstarts, jends):
		wiggle = find_wiggle(start, annotpos[chrom], annotmin[chrom], maxdist)
		if end - start < skipsize:  # junction too small, could be gap
			continue
		if wiggle == maxdist: 
			startstrand = '.'
			gene = ''
		elif start+wiggle in annotpos[chrom]:
			startstrand = '+'
			gene = annotpos[chrom][start+wiggle]
		else:
			startstrand = '-'
			gene = annotmin[chrom][start+wiggle]

		wiggle = find_wiggle(end, annotpos[chrom], annotmin[chrom], maxdist)
		if wiggle == maxdist:
			endstrand = '.'
			endgene = ''
		elif end+wiggle in annotpos[chrom]:
			endstrand = '+'
			endgene = annotpos[chrom][end+wiggle]
		else:
			endstrand = '-'
			endgene = annotmin[chrom][end+wiggle]

		if startstrand == endstrand:  # both +, -, or .
			if startstrand == '.':  # both .
				continue
			consensus = startstrand
		elif startstrand == '.':
			consensus = endstrand
			gene = endgene
		else:
			consensus = startstrand
		strands += [consensus]
		gene_candidates += [gene]
		if len(strands) > skipnum and len(set(strands) - set(['.'])) == 1:
			break
	return strands, gene_candidates

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	dist = {}
	prevline, maxdist = '', 100
	for line in psl:  # compare psl junctions with annotated junctions
		line = line.rstrip().split('\t')
		chrom  = line[13]
		if chrom not in annotpos or chrom not in annotmin:
			continue

		starts = [int(x) for x in line[20].split(',')[:-1]]  # block/exon starts
		ends = [int(x) + y for x,y in zip(line[18].split(',')[:-1], starts)]
		jstarts = ends[:-1]  # junction/intron starts
		jends = starts[1:]
		if len(ends) == 1:
			jstarts = starts
			jends = ends
		else:
			jstarts += [starts[0]]
			jends += [ends[-1]]

		strands, gene_candidates = find_strand(annotpos, annotmin, jstarts, jends, 20, 30)
		pluscount = strands.count('+')
		minuscount = len(strands) - pluscount
		if pluscount == minuscount:  # ambiguous result, do a more sensitive search
			strands, gene_candidates = find_strand(annotpos, annotmin, jstarts, jends, 150, 5, 10)
			pluscount = strands.count('+')
			minuscount = len(strands) - pluscount

		if pluscount == minuscount:
			consensus = '.'
			gene = chrom+':'+str(jstarts[0])
		else:			
			if pluscount > minuscount:
				consensus = '+'
			else:
				consensus = '-'
			gene = gene_candidates[strands.index(consensus)]

		line[8] = consensus
		if not genepred:
			line[9] += '_' + gene
		writer.writerow(line)

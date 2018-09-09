import sys, csv

try:
	psl = open(sys.argv[1])
	mode = sys.argv[2]
	if sys.argv[2] != 'ginormous' and sys.argv[2] != 'comprehensive':
		mode = 'default'  # all typos will report default
	pslout = sys.argv[3]
except:
	print('usage: script.py collapsed.psl [default/comprehensive/ginormous] filtered.psl')
	sys.exit()

def get_junctions_psl(starts, sizes):
	junctions = set()
	if len(starts) != 1:
		for b in range(len(starts)-1):
			junctions.add((starts[b]+sizes[b], starts[b+1]))
		return junctions

def get_junctions_bed12(line):
	junctions = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1):
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_exons_psl(starts, sizes):
	exons = []
	for e in range(len(starts)):
		exons += [(starts[e], starts[e]+sizes[e])]
	return exons

def overlap(coords0, coords1, tol=1):
	coords0, coords1 = sorted([coords0, coords1], key = lambda x: x[0])
	return (coords0[0] < coords1[0] and coords1[0] < coords0[1] - tol)

def contained(coords0, coords1):
	return coords1[0] <= coords0[0] and coords1[1] >= coords0[1]  # complete coverage of coords0 by coords1

def bin_search(query, data):
	""" Query is coordinates. Binary search for the query in data, which is a list of coordinates.
	Finishes when an overlapping value of query and data exists and returns the index in data. """
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

isoforms = {}
allevents = {}
for line in psl:
	line = line.rstrip().split()
	chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])  # psl
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	starts = [int(n) for n in line[20].split(',')[:-1]]
	junctions = get_junctions_psl(starts, sizes)
	exons = get_exons_psl(starts, sizes)
	j = str(sorted(list(junctions)))
	if chrom not in isoforms:
		isoforms[chrom] = {}
		allevents[chrom] = {}
		allevents[chrom]['alljunctions'] = set()
		allevents[chrom]['allexons'] = set()
	if '-' in name[-4:]:
		name = name[:name.rfind('-')]
	if name not in isoforms[chrom]:
		isoforms[chrom][name] = {}
		isoforms[chrom][name]['line'] = [line]
		isoforms[chrom][name]['jset'] = junctions
		isoforms[chrom][name]['exons'] = exons
		allevents[chrom]['alljunctions'].update(junctions)
		allevents[chrom]['allexons'].update(exons)

keepisoforms = []
for chrom in isoforms:
	sys.stderr.write(chrom+'\n')
	alljunctions = sorted(list(allevents[chrom]['alljunctions']), key=lambda x: x[0])
	allexons = sorted(list(allevents[chrom]['allexons']), key=lambda x: x[0])
	for n in isoforms[chrom]:
		if mode == 'ginormous':  # currently includes all
			keepisoforms += isoforms[chrom][n]['line']
			continue
		if len(isoforms[chrom][n]['jset']) >= 1:  # multi-exon isoforms
			if mode == 'comprehensive':  # all spliced subsets allowed
				keepisoforms += isoforms[chrom][n]['line']
				continue
			issubset = set()
			for n_ in isoforms[chrom]:  # compare with all others
				if n == n_:
					continue
				elif isoforms[chrom][n]['jset'] < isoforms[chrom][n_]['jset']:
					issubset.update(isoforms[chrom][n_]['jset'])
					break
			if issubset: 
				ir = False
				for exon in sorted(list(isoforms[chrom][n]['exons'])):
					for junction in sorted(list(issubset)):
						if overlap(exon, junction, 5):
							ir = True
							break
				if ir:  # if there is intron retention, it is not considered a subset isoform
					keepisoforms += isoforms[chrom][n]['line']
			else:  # full set of splice junctions
				keepisoforms += isoforms[chrom][n]['line']
		else:  # single exon isoforms
			exon = list(isoforms[chrom][n]['exons'])[0]
			i = bin_search(exon, alljunctions)
			iscontained = False 
			for j in alljunctions[i-2:i+2]:
				if contained(exon, j):
					iscontained = True
					break
			if iscontained:  # boolean: exon is completely contained within a splice junction
				keepisoforms += isoforms[chrom][n]['line']
				continue

			i = bin_search(exon, allexons)
			iscontained = False  # boolean: exon completely contained in another exon
			for e in allexons[i-2:i+2]:
				if contained(exon, e):
					iscontained = True
					break
			if not iscontained:
				keepisoforms += isoforms[chrom][n]['line']


with open(pslout, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for iso in keepisoforms:
		writer.writerow(iso)

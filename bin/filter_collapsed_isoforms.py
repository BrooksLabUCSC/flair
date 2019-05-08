import sys, csv

try:
	psl = open(sys.argv[1])
	mode = sys.argv[2]
	if sys.argv[2] != 'ginormous' and sys.argv[2] != 'comprehensive':
		mode = 'default'  # all typos will report default
	pslout = sys.argv[3]
	if len(sys.argv) > 4:
		tol = int(sys.argv[4])
	else:
		tol = 10
except:
	print('usage: script.py collapsed.psl (default/comprehensive/ginormous) filtered.psl [tolerance]')
	sys.exit(1)

def get_junctions_psl(starts, sizes):
	junctions = set()
	if len(starts) != 1:
		for b in range(len(starts)-1):
			junctions.add((starts[b]+sizes[b], starts[b+1]))
		return junctions

def get_junction_exon_bed12(line):
	junctions = set()
	exons = set()
	chrstart = int(line[1])
	starts = [int(n) + chrstart for n in line[11].split(',')[:-1]]
	sizes = [int(n) for n in line[10].split(',')[:-1]]
	if len(starts) != 1:
		for b in range(len(starts)-1):
			junctions.add((starts[b]+sizes[b], starts[b+1]))
			exons.add((starts[b], starts[b]+sizes[b]))
		return junctions

def get_exons_psl(starts, sizes):
	exons = []
	for e in range(len(starts)):
		exons += [(starts[e], starts[e]+sizes[e])]
	return exons

def overlap(coords0, coords1, tol=1):
	coords0, coords1 = sorted([coords0, coords1], key = lambda x: x[0])
	return (coords0[0] < coords1[0] and coords1[0] < coords0[1] - tol)

def exon_overlap(coords0, coords1, left=True, tol=1):
	len0, len1 = coords0[1] - coords0[0], coords1[1] - coords1[0]
	if left:
		return coords0[1] == coords1[1] and (len0 - tol) < len1
	return coords0[0] == coords1[0] and (len0 - tol) < len1

def contained(coords0, coords1, tol=0):  # complete coverage of coords0 by coords1
	return coords0[1] > coords1[0] and coords1[0] <= coords0[0]+tol and coords1[1] >= coords0[1]-tol

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
	
isoforms, allevents, jcn_to_name = {}, {}, {}

for line in psl:
	line = line.rstrip().split()
	if sys.argv[1][-3:] != 'psl':  # untested but it's just two lines so should be ok?
		chrom, name = line[0], line[3]
		junctions, exons  = get_junction_exon_bed12(line)
	else:
		chrom, name = line[13], line[9]
		sizes = [int(n) for n in line[18].split(',')[:-1]]
		starts = [int(n) for n in line[20].split(',')[:-1]]
		junctions = get_junctions_psl(starts, sizes)
		exons = get_exons_psl(starts, sizes)
	if chrom not in isoforms:
		isoforms[chrom] = {}
		allevents[chrom] = {}
		allevents[chrom]['alljunctions'] = set()
		allevents[chrom]['allexons'] = set()
		jcn_to_name[chrom] = {}
	if name not in isoforms[chrom]:
		isoforms[chrom][name] = {}
		isoforms[chrom][name]['line'] = [line]
		isoforms[chrom][name]['exons'] = exons
		if junctions:  # multi-exon isoform
			isoforms[chrom][name]['jset'] = junctions
			allevents[chrom]['alljunctions'].update(junctions)
			allevents[chrom]['allexons'].update(exons)
			jname = str(sorted(list(junctions)))
			isoforms[chrom][name]['jname'] = jname
			for j in junctions:
				if j not in jcn_to_name[chrom]:
					jcn_to_name[chrom][j] = set()
				jcn_to_name[chrom][j].add(name)

keepisoforms = []
for chrom in isoforms:
	# sys.stderr.write(chrom+'\n')
	alljunctions = sorted(list(allevents[chrom]['alljunctions']), key=lambda x: x[0])
	allexons = sorted(list(allevents[chrom]['allexons']), key=lambda x: x[0])
	for n in isoforms[chrom]:
		if mode == 'ginormous':  # currently includes all
			keepisoforms += isoforms[chrom][n]['line']
			continue
		if 'jset' in isoforms[chrom][n]:  # multi-exon isoforms
			if mode == 'comprehensive':  # all spliced subsets allowed
				keepisoforms += isoforms[chrom][n]['line']
				continue
			similar_isos = set()
			for j in isoforms[chrom][n]['jset']:
				similar_isos.update(jcn_to_name[chrom][j])

			issubset = [0, 0]  # first exon is a subset, last exon is a subset
			subset_juncs = set()  # if true, will be populated with introns of other isoforms
			exons = sorted(list(isoforms[chrom][n]['exons']),key=lambda x: x[0])
			first_exon, last_exon = exons[0], exons[-1]
			for n_ in similar_isos:  # compare with all similar
				if n == n_:
					continue
				elif isoforms[chrom][n]['jname'][1:-1] in isoforms[chrom][n_]['jname'] and \
				len(isoforms[chrom][n]['jname']) < len(isoforms[chrom][n_]['jname']):
					for exon in list(isoforms[chrom][n_]['exons']):
						if exon_overlap(first_exon, exon, tol=tol):
							issubset[0] = 1
						if exon_overlap(last_exon, exon, left=False, tol=tol):
							issubset[1] = 1
					subset_juncs.update(isoforms[chrom][n_]['jset'])
			if sum(issubset) != 2:
				keepisoforms += isoforms[chrom][n]['line']
		else:  # single exon isoforms
			exon = list(isoforms[chrom][n]['exons'])[0]
			i = bin_search(exon, alljunctions)
			iscontained = False 
			for j in alljunctions[i-2:i+2]:
				if contained(exon, j):
					iscontained = True
					break
			if iscontained:  # boolean: exon is contained within a splice junction
				keepisoforms += isoforms[chrom][n]['line']
				continue
			i = bin_search(exon, allexons)
			iscontained = False  # boolean: exon contained in another exon
			for e in allexons[i-2:i+2]:
				if contained(exon, e, 5):
					iscontained = True
					break
			if not iscontained:
				keepisoforms += isoforms[chrom][n]['line']

with open(pslout, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for iso in keepisoforms:
		writer.writerow(iso)
